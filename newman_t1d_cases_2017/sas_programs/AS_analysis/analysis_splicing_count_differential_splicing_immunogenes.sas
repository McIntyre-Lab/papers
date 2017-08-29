/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Merge in autoimmune gene flags and count fusions for these subsets */

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

data event2gene;
  set splice.splicing_events_annotations;
  keep event_id gene_id;
run;

proc sort data=event2gene nodup;
  by gene_id event_id;
proc sort data=ai nodup;
  by gene_id;
proc sort data=t1d nodup;
  by gene_id;
run;

data splicing_immuno;
  merge event2gene (in=in1) ai (in=in2) t1d (in=in3);
  by gene_id;
  if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
  if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
  if in1 then output;
run;

/* Merge with fusion flags */

data flag_splicing_spec;
   set splicing.flag_splicing_by_gene_dtct_v2;
run;

proc sort data=flag_splicing_spec;
  by event_id;
proc sort data=splicing_immuno;
  by event_id;
run;

data flag_splicing_spec_immuno;
  merge flag_splicing_spec (in=in1) splicing_immuno (in=in2);
  by event_id;
  if in1 and in2;
run;

*check: how many events from autoimmune genes are cell specific?;

proc freq data=flag_splicing_spec_immuno;
   tables flag_cell_specific*flag_immuno_gene;
run;


data check;
  set flag_splicing_spec_immuno;
  where flag_cell_specific=1 and flag_immuno_gene=1;
run; *22 events!;

data check2;
  set flag_splicing_spec_immuno;
  where flag_cell_specific=1 and flag_immunobase_diabetes_gene=1;
run; *3 events!;



/* Count and export */

%macro counts(prefix,outname);

proc freq data=flag_splicing_spec_immuno noprint;
   %if &prefix. = ja %then %do; where flag_junction_annotated=1 and flag_immuno_gene=1; %end;
   %else %if &prefix. = ju %then %do; where flag_junction_annotated=0 and flag_intron_retention=0  and flag_immuno_gene=1; %end;
   %else %if &prefix. = ir %then %do; where flag_intron_retention=1 and flag_immuno_gene=1; %end;
   %else %if &prefix. = es %then %do; where flag_exonskip=1 and flag_immuno_gene=1; %end;
   %else %if &prefix. = ad %then %do; where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_immuno_gene=1; %end;
   %else %if &prefix. = aa %then %do; where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1; %end;
   %else %if &prefix. = ada %then %do; where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immuno_gene=1; %end;
   %else %do; where flag_immuno_gene=1; %end;
   tables flag_cell_specific*flag_alt_spliced*flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4_on*flag_cd8_on*flag_cd19_on*flag_anova_fdr_05*sum_flags / out = &prefix._counts_ai;
run;

proc freq data=flag_splicing_spec_immuno noprint;
   %if &prefix. = ja %then %do; where flag_junction_annotated=1 and flag_immunobase_diabetes_gene=1; %end;
   %else %if &prefix. = ju %then %do; where flag_junction_annotated=0 and flag_intron_retention=0  and flag_immunobase_diabetes_gene=1; %end;
   %else %if &prefix. = ir %then %do; where flag_intron_retention=1 and flag_immunobase_diabetes_gene=1; %end;
   %else %if &prefix. = es %then %do; where flag_exonskip=1 and flag_immunobase_diabetes_gene=1; %end;
   %else %if &prefix. = ad %then %do; where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and flag_immunobase_diabetes_gene=1; %end;
   %else %if &prefix. = aa %then %do; where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immunobase_diabetes_gene=1; %end;
   %else %if &prefix. = ada %then %do; where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and flag_immunobase_diabetes_gene=1; %end;
   %else %do; where flag_immunobase_diabetes_gene=1; %end;

   tables flag_cell_specific*flag_alt_spliced*flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4_on*flag_cd8_on*flag_cd19_on*flag_anova_fdr_05*sum_flags / out = &prefix._counts_t1d;
run;

/* Export data */

data &prefix._counts_ai2;
   set &prefix._counts_ai;
   drop PERCENT;
run;


data &prefix._counts_t1d2;
   set &prefix._counts_t1d;
   drop PERCENT;
run;


/* Export data */

proc export data=&prefix._counts_ai2
     outfile="!PATCON/pipeline_output/differential_splicing_counts/counts_by_event_and_gene_detection_&outname._autoimmune_genes.csv"
     dbms=csv replace;
run;


proc export data=&prefix._counts_t1d2
     outfile="!PATCON/pipeline_output/differential_splicing_counts/counts_by_event_and_gene_detection_&outname._diabetes_genes.csv"
     dbms=csv replace;
run;

%mend;

%counts(all,all_events);
%counts(ja,annotated_junctions);
%counts(ju,unannotated_junctions);
%counts(ir,intron_retention);
%counts(es,exonskipping);
%counts(ad,alt_donor);
%counts(aa,alt_acceptor);
%counts(ada,alt_donor_and_acceptor);

