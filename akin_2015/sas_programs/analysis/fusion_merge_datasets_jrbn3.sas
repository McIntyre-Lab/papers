/* import libraries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';


/* Merge all to table! */

data anova3;
   set sugrue.anova_logapn;
   if HypothesisType=1 then delete;
run;

data transposed_fusion_data;
   set sugrue.fusion_data_flip;
   drop _NAME_;
run;


/* sort all  by fusion_id */
/* first get basic info for sugrue */

data sugrue_fusion_info;
   set sugrue.counts_w_means;
   keep fusion_id region_length flag_control_on flag_treat_on flag_all_on mean_apn_control mean_apn_treat flag_up_down fold_change flag_low_exp_con flag_low_exp_treat flag_low_exp_both;
run;

proc sort data=sugrue_fusion_info nodups;
   by fusion_id;
run;

proc sort data=anova3;
   by fusion_id;
run;

proc sort data=transposed_fusion_data;
   by fusion_id;
run;

proc sort data=sugrue.flag_resids_logapn;
   by fusion_id;
run;

data anova_resid oops1 oops2;
  merge anova3 (in=in1) sugrue.flag_resids_logapn (in=in2);
  by fusion_id;
  if in1 and in2 then output anova_resid;
  else if in1 then output oops1;
  else output oops2;
run;

data anova_resid_counts oops;
  merge transposed_fusion_data (in=in1) anova_resid (in=in2);
  by fusion_id;
  if in1 and in2 then output anova_resid_counts;
  else if in1 then output anova_resid_counts; *missing anova results will be set to missing for non-analyzed fusions;
  else output oops;
run;

data results_by_fusion oops1 oops2;
  merge sugrue_fusion_info (in=in1) anova_resid_counts (in=in2);
  by fusion_id;
  if in1 and in2 then output results_by_fusion;
  else if in1 then output oops1;
  else output oops2;
run;

/* make permenant */

data sugrue.results_by_fusion;
   set results_by_fusion;
   if  ProbF=. then flag_p05=.;
   else if ProbF < 0.05 then flag_p05=1;
   else flag_p05=0;
run;


