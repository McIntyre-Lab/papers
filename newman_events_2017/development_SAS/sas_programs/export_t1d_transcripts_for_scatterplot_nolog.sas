ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For transcripts in both lists, export variance and make a scatterplot */

data xs2keep;
   set event.hg19_flag_xscripts_w_unique;
   if perc_features_dtct_cd4 >= 0.75
   or perc_features_dtct_cd8 >= 0.75
   or perc_features_dtct_cd19 >= 0.75;
   keep transcript_id;
run;

data all_var;
  set eventloc.hg19_variance_all_xs_nolog;
  keep transcript_id var_cd19 var_cd4 var_cd8;
  rename var_cd19=var_cd19_all var_cd4=var_cd4_all var_cd8=var_cd8_all;
run;

data subset_var;
   set eventloc.hg19_variance_filtered_xs_nolog;
  keep transcript_id var_cd19 var_cd4 var_cd8;
  rename var_cd19=var_cd19_subset var_cd4=var_cd4_subset var_cd8=var_cd8_subset;
run;

proc sort data=xs2keep;
   by transcript_id;
proc sort data=all_var;
   by transcript_id;
proc sort data=subset_var;
   by transcript_id;
run;


data var_all_v_subset;
  merge xs2keep (in=in1) all_var (in=in2) subset_var (in=in3);
  by transcript_id;
  if in1 and in2;
run;

proc export data=var_all_v_subset
     outfile="!MCLAB/event_analysis/analysis_output/hg19_xscript_variances_common_nolog.csv"
     dbms=csv replace;
run;

