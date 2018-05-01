ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';



/* Export data for making plots to look at transcript binning */

proc export data=event.hg19_flag_xscripts_w_unique
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_t1d_transcripts_by_perc_unique_detected_nomulti.csv" dbms=csv replace;
run;

