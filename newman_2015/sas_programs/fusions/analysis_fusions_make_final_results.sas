/* import libraries */
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/*** Merge datasets ***/

/* need:
residual normality
transposed data
means 
FDR
fusion_info
*/

/* sort all  by fusion_id */
/* first get basic info for sugrue */

data fusion_info;
   set fus.unique_info_fusions_si;
   keep fusion_id gene_id num_exons flag_multigene;
run;

proc sort data=fusion_info nodups;
   by fusion_id;
run;

proc sort data=con.results_by_fusion_w_fdr;
   by fusion_id;
run;

proc sort data=Con.fusions_resid_normtest;
   by fusion_id;


proc sort data=con.fusion_means_by_celltype;
  by fusion_id;
run;

proc sort data=con.fusions_on_gt_apn0;
  by fusion_id;
run;


data fusion_w_flags;
   merge fusion_info (in=in1) con.fusions_on_gt_apn0 (in=in2);
   by fusion_id;
   if flag_CD19_on=. then flag_CD19_on=0;
   if flag_CD4_on=. then flag_CD4_on=0;
   if flag_CD8_on=. then flag_CD8_on=0;
   if flag_fusion_on0=. then flag_fusion_on0=0;
   if flag_fusion_all_on0=. then flag_fusion_all_on0=0;
run;

data fusion_w_means;
   merge fusion_w_flags (in=in1)  con.fusion_means_by_celltype (in=in2);
   by fusion_id;
   if mean_logq3q3apn_cd19=. then mean_logq3q3apn_cd19=0;
   if mean_logq3q3apn_cd4=. then mean_logq3q3apn_cd4=0;
   if mean_logq3q3apn_cd8=. then mean_logq3q3apn_cd8=0;
run;

data fusion_w_normtest;
   merge fusion_w_means (in=in1) Con.fusions_resid_normtest (in=in2);
   by fusion_id;
   if pnorm=. then flag_fail_norm=.;
   else if pnorm lt 0.05 then flag_fail_norm=1;
   else flag_fail_norm=0;
run;

data fusion_w_results;
   merge fusion_w_normtest (in=in1) con.results_by_fusion_w_fdr (in=in2);
   by fusion_id;
run;

/* Make permenant */

data con.results_by_fusion_final;
   set fusion_w_results;
run;

/* Export results file */

proc export data=con.results_by_fusion_final
    outfile='/home/jrbnewman/concannon/generated_files/results_by_fusion.csv'
    dbms=csv
    replace;
    putnames=yes;
run;
