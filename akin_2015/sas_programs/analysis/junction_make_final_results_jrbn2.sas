/* JUNCTIONS */

/* Because the junctions datasets are large I will be working on this locally. I will save the final datasets to the share */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';


/* Merge all to table! */

data anova3;
   set mysas.junc_anova_&sysparm.;
   if HypothesisType=1 then delete;
   keep fusion_id DF FValue ProbF;
run;

/* sort all  by fusion_id */

data sugrue_junction_info;
   set mysas.junc_counts_w_annot_&sysparm.;
   drop _NAME_;
run;

proc sort data=sugrue_junction_info nodups;
   by fusion_id;
run;

proc sort data=anova3;
   by fusion_id;
run;

proc sort data=mysas.junc_resid_flags_&sysparm.;
   by fusion_id;
run;

data anova_resid oops1 oops2;
  merge anova3 (in=in1) mysas.junc_resid_flags_&sysparm. (in=in2);
  by fusion_id;
  if in1 and in2 then output anova_resid;
  else if in1 then output oops1;
  else output oops2;
run;

data results_by_junction oops2;
  merge sugrue_junction_info (in=in1) anova_resid (in=in2);
  by fusion_id;
  if in1 and in2 then output results_by_junction;
  else if in1 then output results_by_junction;
  else output oops2;
run;

/* make permenant */

data mysas.results_by_junction_&sysparm.;
   set results_by_junction;
   length flag_up_down_exp $1.;
   if ProbF=. then flag_p05=.;
   else if ProbF < 0.05 then flag_p05=1;
   else flag_p05=0;
   if flag_p05=1 then flag_up_down_exp=flag_up_down;
   else flag_up_down_exp='.';
run;

