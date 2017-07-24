ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For transcripts in both lists, export variance and make a scatterplot */

data all_var;
  set eventloc.hg19_cv_all_xs;
  keep transcript_id cv_cd19 cv_cd4 cv_cd8;
  rename cv_cd19=cv_cd19_all cv_cd4=cv_cd4_all cv_cd8=cv_cd8_all;
run;

data subset_var;
   set eventloc.hg19_cv_filtered_xs;
  keep transcript_id cv_cd19 cv_cd4 cv_cd8;
  rename cv_cd19=cv_cd19_subset cv_cd4=cv_cd4_subset cv_cd8=cv_cd8_subset;
run;

proc sort data=all_var;
   by transcript_id;
proc sort data=subset_var;
   by transcript_id;
run;


data var_all_v_subset;
  merge all_var (in=in1) subset_var (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc export data=var_all_v_subset
     outfile="!MCLAB/event_analysis/analysis_output/hg19_xscript_cv_common.csv"
     dbms=csv replace;
run;

