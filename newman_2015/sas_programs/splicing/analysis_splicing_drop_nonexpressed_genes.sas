/* Drop splicing events from genes without at least one exon expressed */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Get list of genes not expressed */

*flag fusions as being on or not on - going to sum these and any genes with a zero will be dropped;

data fusion_exp;
  set con.fusions_on_gt_apn0;
  if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=0 then flag_fusion_on=0;
  else flag_fusion_on=1;
  keep fusion_id flag_fusion_on;
run;


*Since not dealing with multigene splicing events, need to summarize fusions info to individual genes, and not multigenes;
data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   keep fusion_id gene_id;
run;

proc sort data=fus2gene;
   by fusion_id;
proc sort data=fusion_exp;
   by fusion_id;
run;

data fus2gene_w_exp oops;
   merge fus2gene (in=in1) fusion_exp (in=in2);
   by fusion_id;
   if in1 and in2 then output fus2gene_w_exp;
   else output oops;
run;

*Summarize expression flag on gene;

proc sort data=fus2gene_w_exp;
   by gene_id;
run;

proc means data=fus2gene_w_exp noprint;
   by gene_id;
   var flag_fusion_on;
   output out=gene_exp sum=num_exons_exp;
run;

data genes_off;
   set gene_exp;
   if num_exons_exp=0;
run;

data splicing.genes_off;
set genes_off;
run;

/* Drop splicing events */

data splicing_results;
   set splicing.splicing_results_w_annot_fdr;
   keep event_id flag_cd19_on flag_cd4_on flag_cd8_on flag_anova_fdr_05 
        gene_id flag_junction_annotated flag_intron_retention flag_exonskip flag_alt_donor flag_alt_acceptor;
   run;

data genes_off;
   set splicing.genes_off;
   keep gene_id;
run;

proc sort data=splicing_results;
   by gene_id;
proc sort data=genes_off;
   by gene_id;
run;

data splicing_gene_on splicing_gene_off;
   merge genes_off (in=in1) splicing_results;
   by gene_id;
   if in1 then output splicing_gene_off;
   else output splicing_gene_on;
run;

/* Make clean dataset permenant -- we will recalc everything off of this */

data splicing.splicing_results_clean;
   set splicing_gene_on;
run;

