libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

data xs_list;
  set event.polyester_xs_list_60genes;
  gene_id2=gene_id+0;
  drop gene_id;
  rename gene_id2=gene_id;
run;

proc import datafile="/home/jrbnewman/McLab/useful_mouse_data/mm10/gff/mm10_refseq_entrezID_symbol.txt"
     out=sym dbms=tab replace;
     guessingrows=max; getnames=no;
run;

data sym2;
  set sym;
  rename VAR1=gene_id VAR2=symbol;
run;

proc sort data=xs_list;
  by gene_id;
proc sort data=sym2;
  by gene_id;
run;

data xs_list_w_sym;
  merge xs_list (in=in1) sym2 (in=in2);
  by gene_id;
  if in1;
run;

proc sort data=xs_list_w_sym;
  by gene_id transcript_id;
proc freq data=xs_list_w_sym noprint;
   tables gene_id / out=xs_count;
proc sort data=xs_count;
   by descending count;
proc print data=xs_count (obs=1);
run; *10 x??;



data cat_xs;
   array xs[10] $15.;
   retain xs1-xs10;
   set xs_list_w_sym;
   by gene_id;
   if first.gene_id then do;
      call missing(of xs1-xs10);
      records = 0;
   end;
   records + 1;
   xs[records]=transcript_id;
   if last.gene_id then output;
run;


data cat_xs2;
  set cat_xs;
  length transcript_id2 $1500.;
  transcript_id2=catx(",", OF xs1-xs10);
  drop xs1-xs10 transcript_id records;
  rename transcript_id2=transcript_id;
run;

proc export data=cat_xs2 outfile="!MCLAB/event_analysis/analysis_output/59_genes_and_transcripts_list.txt"
    dbms=tab replace;
run;



libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

data xs_list;
  set event.polyester_xs_list_10k;
run;

data xs2gene;
  set event.feature2xs2gene;
  gene_id2=gene_id+0;
  keep transcript_id gene_id2;
  rename gene_id2=gene_id;
run;

proc sort data=xs_list;
  by transcript_id;
proc sort data=xs2gene nodup;
  by transcript_id gene_id;
run;

data xs_list2;
  merge xs2gene (in=in1) xs_list (in=in2);
  by transcript_id;
  if in1 and in2;
run;


proc import datafile="/home/jrbnewman/McLab/useful_mouse_data/mm10/gff/mm10_refseq_entrezID_symbol.txt"
     out=sym dbms=tab replace;
     guessingrows=max; getnames=no;
run;

data sym2;
  set sym;
  rename VAR1=gene_id VAR2=symbol;
run;

proc sort data=xs_list2;
  by gene_id;
proc sort data=sym2;
  by gene_id;
run;

data xs_list_w_sym;
  merge xs_list2 (in=in1) sym2 (in=in2);
  by gene_id;
  if in1;
run;

proc export data=xs_list_w_sym outfile="!MCLAB/event_analysis/analysis_output/10000_transcripts_list.txt"
    dbms=tab replace;
run;



