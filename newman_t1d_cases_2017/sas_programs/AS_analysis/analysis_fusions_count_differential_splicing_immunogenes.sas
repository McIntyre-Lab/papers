
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


data fus2gene;
  set fus.hg19_aceview_fusions_si_info;
  keep fusion_id gene_id;
run;


proc sort data=fus2gene nodup;
  by gene_id fusion_id;
proc sort data=ai nodup;
  by gene_id;
proc sort data=t1d nodup;
  by gene_id;
run;

data fusion_immuno;
  merge fus2gene (in=in1) ai (in=in2) t1d (in=in3);
  by gene_id;
  if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
  if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
run;

proc sort data=fusion_immuno nodup;
   by fusion_id;
proc means data=fusion_immuno noprint;
   by fusion_id;
   var flag_immuno_gene flag_immunobase_diabetes_gene;
   output out=fusion_immuno2 max=;
run;


/* Merge with fusion flags */

data flag_fusion_spec;
   set con.flag_fusions_by_gene_detection2;
run;

proc sort data=flag_fusion_spec;
  by fusion_id;
proc sort data=fusion_immuno2;
  by fusion_id;
run;

data flag_fusion_spec_immuno;
  merge flag_fusion_spec (in=in1) fusion_immuno2 (in=in2);
  by fusion_id;
  if in1 and in2;
  drop _TYPE_ _FREQ_;
run;

/* Count and export */

proc freq data=flag_fusion_spec_immuno noprint;
   where flag_immuno_gene=1;
   tables flag_cell_specific*flag_alt_spliced*
          flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4_on*flag_cd8_on*flag_cd19_on*
          flag_cd4cd8_fdr05*flag_cd4cd19_fdr05*flag_cd8cd19_fdr05*sum_flags / out = fusion_counts_ai;
run;

proc freq data=flag_fusion_spec_immuno noprint;
   where flag_immunobase_diabetes_gene=1;
   tables flag_cell_specific*flag_alt_spliced*
          flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on*
          flag_cd4_on*flag_cd8_on*flag_cd19_on*
          flag_cd4cd8_fdr05*flag_cd4cd19_fdr05*flag_cd8cd19_fdr05*sum_flags / out = fusion_counts_t1d;
run;

/* Make dataset permenant */


/* Export data */

data fusion_counts_ai2;
   set fusion_counts_ai;
   drop PERCENT;
run;

proc export data=fusion_counts_ai2
     outfile="!PATCON/pipeline_output/differential_splicing_counts/counts_by_fusion_and_gene_detection_autoimmune_genes.csv"
     dbms=csv replace;
run;


data fusion_counts_t1d2;
   set fusion_counts_t1d;
   drop PERCENT;
run;

proc export data=fusion_counts_t1d2
     outfile="!PATCON/pipeline_output/differential_splicing_counts/counts_by_fusion_and_gene_detection_diabetes_genes.csv"
     dbms=csv replace;
run;


