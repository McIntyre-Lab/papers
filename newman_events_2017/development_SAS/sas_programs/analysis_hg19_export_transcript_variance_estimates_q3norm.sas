ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Export variance estimates for plotting the distribution */


data xs2keep;
   set event.hg19_flag_xscripts_w_unique;
   if perc_features_dtct_cd4 ge 0.75
   or perc_features_dtct_cd8 ge 0.75
   or perc_features_dtct_cd19 ge 0.75;
   keep transcript_id;
run;

* all transcripts ;

data flag_on_all;
   set eventloc.flag_xscript_tpm_gt0_all;
run;

data var_all;
  set eventloc.hg19_variance_all_xs_q3norm;
  keep transcript_id var_cd19 var_cd4 var_cd8;
run;

proc sort data=xs2keep;
   by transcript_id;
proc sort data=var_all;
   by transcript_id;
proc sort data=flag_on_all;
   by transcript_id;
run;

data xs_all_data;
  merge var_all (in=in1) xs2keep (in=in2) flag_on_all (in=in3);
  by transcript_id;
  if in2 then flag_filter=1; else flag_filter=0;
  if in1;
run;

* filtered transcripts ;
data flag_on_filt;
   set eventloc.flag_xscript_tpm_gt0_filt;
run;

data var_filtered;
   set eventloc.hg19_variance_filtered_xs_q3norm;
  keep transcript_id var_cd19 var_cd4 var_cd8;
  rename var_cd4=var_cd4_filter var_cd8=var_cd8_filter var_cd19=var_cd19_filter;
run;

proc sort data=xs2keep;
   by transcript_id;
proc sort data=var_filtered;
   by transcript_id;
proc sort data=flag_on_filt;
   by transcript_id;
run;

data xs_filt_data;
   merge xs2keep (in=in1) var_filtered (in=in2) flag_on_filt (in=in3);
   by transcript_id;
   if in1 and in2;
run;

/* Export for plots */

proc export data=xs_all_data
     outfile="!MCLAB/event_analysis/analysis_output/hg19_aceview_all_transcripts_variance_q3norm.csv"
     dbms=csv replace;
run;

proc export data=xs_filt_data
     outfile="!MCLAB/event_analysis/analysis_output/hg19_aceview_75perc_apn5_transcripts_variance_q3norm.csv"
     dbms=csv replace;
run;

