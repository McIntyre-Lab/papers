ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';

/* For transcripts with at least 1 event detected and 1 unique event detected, I want to plot the distribution of % unique detected for each bin: 1-25%, 25-50%, 50-75%, 75-99%, 100%.

So, I am going to export the necessary data, then make the plot in python */

data xscripts_for_export;
  set event.xscripts_w_unique_by_bin_total;
  if perc_unique_features_dtct=. then delete;
  if perc_unique_features_dtct=0 then delete;
  keep transcript_id perc_features_dtct perc_unique_features_dtct;
run;


proc export data=xscripts_for_export outfile="!MCLAB/event_analysis/xscripts_w_uniq_dtct_total_dtct.csv"
   dbms=csv replace;
run;
