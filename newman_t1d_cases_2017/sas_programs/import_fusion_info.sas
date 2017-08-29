
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* import fusion files into SAS */ 

%macro import_fusions(type);

proc import datafile="/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/generated_files/hg19_aceview_fusions_&type..bed" out=fusions_&type._bed dbms=tab replace;
   getnames=no; guessingrows=349000;
   run;

proc import datafile="/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/generated_files/hg19_aceview_fusions_&type..tsv" out=fusions_&type._info dbms=tab replace;
   getnames=yes; guessingrows=680000;
   run;

/* Make permenant */

data fus.hg19_aceview_fusions_&type._bed;
  set fusions_&type._bed;
  rename VAR1=chr
         VAR2=fusion_start
         VAR3=fusion_stop
         VAR4=fusion_id
         VAR5=score
         VAR6=strand;
run;


data fus.hg19_aceview_fusions_&type._info;
  set fusions_&type._info;
  rename primary_FBgn=gene_id;
run;

/* Merge BED and INFO files */

proc sort data=fus.hg19_aceview_fusions_&type._bed;
   by fusion_id;
proc sort data=fus.hg19_aceview_fusions_&type._info;
   by fusion_id;
run;

data fus.exons2fusions_&type. oops1 oops2;
   merge fus.hg19_aceview_fusions_&type._bed (in=in1) fus.hg19_aceview_fusions_&type._info (in=in2);
   by fusion_id;
   if in1 and in2 then output fus.exons2fusions_&type.;
   else if in1 then output oops1;
   else output oops2;
run;

/* Concatenate gene_ids for each fusion */
/* Need BED contents, flag_multigene, num_genes, gene_id, num_exons */

* Exon counts;

data fusion2exon;
   set fus.exons2fusions_&type.;
   keep fusion_id exon_id;
run;

proc sort data=fusion2exon;
   by fusion_id;
run;

proc freq data=fusion2exon noprint;
   table fusion_id / out=num_exons;
run;

* Cat gene_ids;

data fusion2gene;
   set fus.exons2fusions_&type.;
   keep fusion_id gene_id;
run;

proc sort data=fusion2gene nodup;
   by fusion_id gene_id;
run;

*6 genes SI;
*3 genes SD;

data cat_genes; 
  array gene[6] $ 36.;
  retain gene1-gene6;
  set fusion2gene;
  by fusion_id;
  if first.fusion_id then do;
     call missing(of gene1-gene6);
     records = 0;
  end;
  records + 1;
  gene[records]=gene_id;
  if last.fusion_id then output;
run;

  *clean up the output file;
data cat_genes2;
  set cat_genes;
  length gene_cat $231.;
  rename records= num_genes;
         gene_cat= catx("|", OF gene1-gene6);
  drop gene1-gene6 gene_id;
  run;

/* Merge with BED dataset */

proc sort data=cat_genes2;
   by fusion_id;
proc sort data=num_exons;
   by fusion_id;
proc sort data=fus.hg19_aceview_fusions_&type._bed;
   by fusion_id;
run;


data unique_fusion_info oops1 oops2 oops3;
  merge fus.hg19_aceview_fusions_&type._bed (in=in1) num_exons (in=in2) cat_genes2 (in=in3);
  by fusion_id;
  if in1 and in2 and in3 then output unique_fusion_info;
  else do;
      if in1 then output oops1;
      if in2 then output oops2;
      if in2 then output oops3;
      end;
run;

/* Make permenant */

data fus.unique_info_fusions_&type.;
   retain chr fusion_start fusion_stop fusion_id strand flag_multigene num_genes gene_cat COUNT;
   set unique_fusion_info;
   if num_genes gt 1 then flag_multigene=1;
   else flag_multigene=0;
   drop score PERCENT;
   rename gene_cat=gene_id;
   rename strand=fusion_strand; 
   rename COUNT=num_exons;
run;

%mend;

%import_fusions(si);
%import_fusions(sd);


