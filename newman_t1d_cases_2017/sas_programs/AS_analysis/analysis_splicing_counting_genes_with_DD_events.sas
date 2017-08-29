
/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* I want to now count the number of genes detected in each cell type pair that have at least 1 fusion
   differentially detected. This will be evidence of alternative splice (gene must be detected in both cell types!)
*/


data flag_splicing_spec ;
   set splicing.flag_splicing_by_gene_dtct;
run;

data event2gene;
  set splice.splicing_events_annotations;
  keep event_id gene_id;
run;


proc sort data=flag_splicing_spec;
   by event_id;
proc sort data=event2gene nodup;
   by event_id gene_id;
run;

data splicing_spec_w_gene;
  merge flag_splicing_spec (in=in1) event2gene (in=in2);
  by event_id;
  if in1 and in2;
run;

data check;
  set splicing_spec_w_gene;
  if sum(flag_CD19_gene_on,flag_cd4_gene_on,flag_cd8_gene_on) gt 1 and
     sum(flag_CD19_on,flag_cd4_on,flag_cd8_on) gt 0;
run;

data genes;
  set check;
  keep gene_id;
run;

proc sort data=genes nodup;
  by gene_id;
run;

209383 events
 18333 genes
