ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Export variance estimates for plotting the distribution */

data xs2keep;
   set eventloc.hg19_cv_filtered_xs;
   keep transcript_id;
run;

* all transcripts ;

data var_all;
  set eventloc.hg19_cv_all_xs;
  keep transcript_id cv_cd19 cv_cd4 cv_cd8;
run;

proc sort data=xs2keep;
   by transcript_id;
proc sort data=var_all;
   by transcript_id;
run;

data xs_all_data;
  merge var_all (in=in1) xs2keep (in=in2);
  by transcript_id;
  if in2 then flag_filter=1; else flag_filter=0;
  if in1;
run;

* filtered transcripts ;

data var_filtered;
   set eventloc.hg19_cv_filtered_xs;
  keep transcript_id cv_cd19 cv_cd4 cv_cd8;
  rename cv_cd4=cv_cd4_filter cv_cd8=cv_cd8_filter cv_cd19=cv_cd19_filter;
run;

/* Export for plots */

proc export data=xs_all_data
     outfile="!MCLAB/event_analysis/analysis_output/hg19_aceview_all_transcripts_cv.csv"
     dbms=csv replace;
run;

proc export data=var_filtered
     outfile="!MCLAB/event_analysis/analysis_output/hg19_aceview_75perc_apn5_transcripts_cv.csv"
     dbms=csv replace;
run;

