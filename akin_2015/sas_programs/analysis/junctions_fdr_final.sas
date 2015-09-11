/* JUNCTIONS */

/* Because the junctions datasets are large I will be working on this locally. I will save the final datasets to the share */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';

data mysas.results_by_junction;
   set mysas.results_by_junction_1 mysas.results_by_junction_2 mysas.results_by_junction_3 mysas.results_by_junction_4 mysas.results_by_junction_5 
mysas.results_by_junction_6 mysas.results_by_junction_7 mysas.results_by_junction_8 mysas.results_by_junction_9 mysas.results_by_junction_10 
mysas.results_by_junction_11 mysas.results_by_junction_12 mysas.results_by_junction_13 mysas.results_by_junction_14 mysas.results_by_junction_15 
mysas.results_by_junction_16 mysas.results_by_junction_17 mysas.results_by_junction_18 mysas.results_by_junction_19 mysas.results_by_junction_20 
mysas.results_by_junction_21 mysas.results_by_junction_22 mysas.results_by_junction_23 mysas.results_by_junction_24 mysas.results_by_junction_25 ;
if flag_all_on=0 then do;
DF=.;
FValue=.;
ProbF=.;
pnorm=.;
flag_fail_norm=.;
flag_up_down_exp=.;
flag_p05=.;
end;
run;


data jnc_for_fdr;
   set mysas.results_by_junction;
   if flag_all_on=1 then output;
   keep fusion_id ProbF;
run;


proc multtest inpvalues(ProbF)=jnc_for_fdr fdr
 out=results_jnct_fdr noprint;
run;
quit;


proc sort data=results_jnct_fdr;
 by fusion_id;
run;

proc sort data=mysas.results_by_junction;
 by fusion_id;
run;

data results_by_jnct_fdr;
  merge mysas.results_by_junction results_jnct_fdr;
  by fusion_id;
run;

proc sort data=results_by_jnct_fdr;
   by fusion_id;
run;

data junc_annot;
   set mysas.splicing_events_annotations;
   rename event_id=fusion_id;
run;

proc sort data=junc_annot;
   by fusion_id;
run;

data mysas.results_by_jnct_fdr;
   merge junc_annot (in=in1) results_by_jnct_fdr;
   by fusion_id;
   if in1;
run;


/* subset data and save to share */

data mysas.exp_jnc_results;
   set mysas.results_by_jnct_fdr;
   if num_samples_exp ge 1 then output;
   exp_indicator=flag_control_on+flag_treat_on;
   drop mean_log_apn_con mean_apn_con mean_log_apn_treat mean_apn_treat flag_goi;
   rename flag_control_on=flag_con_on_50pc;
   rename flag_treat_on=flag_treat_on_50pc;
run;

proc sort data=mysas.exp_jnc_results;
    by flag_p05 gene_id;
run;


data mysas.exp_jnc_analyzed;
   set mysas.exp_jnc_results;
   if flag_all_on = 1 then output;
run;

/* 12506710 splicing events total, 63551 "expressed */

data mysas.goi_jnc_results;
   set mysas.results_by_jnct_fdr;
   if flag_goi=1 then output;
   exp_indicator=flag_control_on+flag_treat_on;
   drop mean_log_apn_con mean_apn_con mean_log_apn_treat mean_apn_treat flag_goi;
   rename flag_control_on=flag_con_on_50pc;
   rename flag_treat_on=flag_treat_on_50pc;
run;

proc sort data=mysas.goi_jnc_results;
    by flag_p05 gene_id;
run;


/* 12506710 splicing events total, 25115 in gene of interest */

data mysas.present_jnc_results;
   set mysas.results_by_jnct_fdr;
   if num_samples_present ge 1 then output;
   exp_indicator=flag_control_on+flag_treat_on;
   drop mean_log_apn_con mean_apn_con mean_log_apn_treat mean_apn_treat flag_goi;
run;


proc sort data=mysas.present_jnc_results;
    by gene_id;
run;

/* 12506710 splicing events total, 63551 "present" - all the same as expressed (not a code bug!) */

data mysas.unannot_jnc_results;
   set mysas.results_by_jnct_fdr;
   if num_samples_exp=0 then delete;
   if flag_junction_annotated=1 then delete;
   exp_indicator=flag_control_on+flag_treat_on;
   drop mean_log_apn_con mean_apn_con mean_log_apn_treat mean_apn_treat flag_goi;
   rename flag_control_on=flag_con_on_50pc;
   rename flag_treat_on=flag_treat_on_50pc;
run;

/* 12506710 splicing events total, 4240 are unannotated and "present" */

proc sort data=mysas.unannot_jnc_results;
    by gene_id;
run;

proc export data=mysas.exp_jnc_results
	outfile='/scratch/lfs/sugrue/sas_analysis/exp_jnc_results2.csv'
	dbms=csv replace;
	run;


proc export data=mysas.exp_jnc_analyzed
	outfile='/scratch/lfs/sugrue/sas_analysis/exp_jnc_analyzed.csv'
	dbms=csv replace;
	run;



proc export data=mysas.unannot_jnc_results
	outfile='/scratch/lfs/sugrue/sas_analysis/unannot_jnc_results2.csv'
	dbms=csv replace;
	run;


proc export data=mysas.present_jnc_results
	outfile='/scratch/lfs/sugrue/sas_analysis/present_jnc_results2.csv'
	dbms=csv replace;
	run;


proc export data=mysas.goi_jnc_results
	outfile='/scratch/lfs/sugrue/sas_analysis/goi_jnc_results2.csv'
	dbms=csv replace;
	run;


