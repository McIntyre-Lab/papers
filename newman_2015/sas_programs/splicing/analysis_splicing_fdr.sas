/*** CALCULATE FDR ***/

/* import libraries */

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data';

/* Cat datasets */

data all_results;
  set mysas.splicing_anova_results_1 mysas.splicing_anova_results_2 mysas.splicing_anova_results_3
  mysas.splicing_anova_results_4 mysas.splicing_anova_results_5 mysas.splicing_anova_results_6 mysas.splicing_anova_results_7
  mysas.splicing_anova_results_8 mysas.splicing_anova_results_9 mysas.splicing_anova_results_10
  mysas.splicing_anova_results_11 mysas.splicing_anova_results_12 mysas.splicing_anova_results_13 
mysas.splicing_anova_results_14 mysas.splicing_anova_results_15 mysas.splicing_anova_results_16 
mysas.splicing_anova_results_17 mysas.splicing_anova_results_18 mysas.splicing_anova_results_19 
mysas.splicing_anova_results_20 mysas.splicing_anova_results_21 mysas.splicing_anova_results_22 
mysas.splicing_anova_results_23 mysas.splicing_anova_results_24 mysas.splicing_anova_results_25;
run;

/* calculate FDR */

proc multtest inpvalues(CD4_CD8_P)=all_results fdr
 out=cd4_8_fdr noprint;
run;
quit;

data cd4_8_fdr2;
  set cd4_8_fdr;
  if fdr_p lt 0.05 then flag_cd4cd8_fdr05=1;
  else flag_cd4cd8_fdr05=0;
  rename fdr_p=CD4_CD8_FDR;
run; 


proc multtest inpvalues(CD4_CD19_P)=cd4_8_fdr2 fdr
 out=cd4_19_fdr noprint;
run;
quit;

data cd4_19_fdr2;
  set cd4_19_fdr;
  if fdr_p lt 0.05 then flag_cd4cd19_fdr05=1;
  else flag_cd4cd19_fdr05=0;
  rename fdr_p=CD4_CD19_FDR;
run; 


proc multtest inpvalues(CD8_CD19_P)=cd4_19_fdr2 fdr
 out=cd8_19_fdr noprint;
run;
quit;

data cd8_19_fdr2;
  set cd8_19_fdr;
  if fdr_p lt 0.05 then flag_cd8cd19_fdr05=1;
  else flag_cd8cd19_fdr05=0;
  rename fdr_p=CD8_CD19_FDR;
run; 

/* Add means and flags back in */

data merged_means;
   set mysas.splicing_means_by_celltype_1 mysas.splicing_means_by_celltype_2 mysas.splicing_means_by_celltype_3
mysas.splicing_means_by_celltype_4 mysas.splicing_means_by_celltype_5 mysas.splicing_means_by_celltype_6
mysas.splicing_means_by_celltype_7 mysas.splicing_means_by_celltype_8 mysas.splicing_means_by_celltype_9
mysas.splicing_means_by_celltype_10 mysas.splicing_means_by_celltype_11 mysas.splicing_means_by_celltype_12 
mysas.splicing_means_by_celltype_13 mysas.splicing_means_by_celltype_14 mysas.splicing_means_by_celltype_15 
mysas.splicing_means_by_celltype_16 mysas.splicing_means_by_celltype_17 mysas.splicing_means_by_celltype_18 
mysas.splicing_means_by_celltype_19 mysas.splicing_means_by_celltype_20 mysas.splicing_means_by_celltype_21 
mysas.splicing_means_by_celltype_22 mysas.splicing_means_by_celltype_23 mysas.splicing_means_by_celltype_24 
mysas.splicing_means_by_celltype_25;

   keep event_id mean_depth_CD19 mean_depth_CD4 mean_depth_CD8;
   run;


data merged_flags;
  set mysas.counts_by_splicing_w_flags_1 mysas.counts_by_splicing_w_flags_2 mysas.counts_by_splicing_w_flags_3
mysas.counts_by_splicing_w_flags_4 mysas.counts_by_splicing_w_flags_5 mysas.counts_by_splicing_w_flags_6
mysas.counts_by_splicing_w_flags_7 mysas.counts_by_splicing_w_flags_8 mysas.counts_by_splicing_w_flags_9
mysas.counts_by_splicing_w_flags_10 mysas.counts_by_splicing_w_flags_11 mysas.counts_by_splicing_w_flags_12
mysas.counts_by_splicing_w_flags_13 mysas.counts_by_splicing_w_flags_14 mysas.counts_by_splicing_w_flags_15 
mysas.counts_by_splicing_w_flags_16 mysas.counts_by_splicing_w_flags_17 mysas.counts_by_splicing_w_flags_18 
mysas.counts_by_splicing_w_flags_19 mysas.counts_by_splicing_w_flags_20 mysas.counts_by_splicing_w_flags_21 
mysas.counts_by_splicing_w_flags_22 mysas.counts_by_splicing_w_flags_23 mysas.counts_by_splicing_w_flags_24 
mysas.counts_by_splicing_w_flags_25 ;
  keep event_id flag_CD19_on flag_CD4_on flag_CD8_on flag_all_on count_CD19 count_CD4 count_CD8;
  run;

proc sort data=merged_means nodup;
   by event_id;
run;

proc sort data=merged_flags nodup;
   by event_id;
run;


proc sort data=cd8_19_fdr2;
  by event_id;
run;


data results_by_splicing_w_fdr;
   merge merged_means merged_flags cd8_19_fdr2;
   by event_id;
run;


/* Make permenant */

data mysas.results_by_splicing_w_fdr;
   set results_by_splicing_w_fdr;
run;


