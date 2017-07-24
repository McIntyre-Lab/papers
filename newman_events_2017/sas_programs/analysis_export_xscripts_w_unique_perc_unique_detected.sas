ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For transcripts with at least one unique event, I want to plot the number of total unique events
   against the proportion of unique events detected */

data xs_uniq;
   set event.xscripts_w_unique_by_bin;
   where flag_xscript_has_unique=1;
   keep transcript_id num_unique_features perc_unique_features_dtct;
run;

proc export data=xs_uniq
     outfile="!MCLAB/event_analysis/analysis_output/xscripts_uniq_feat_by_perc_dtct.csv"
     dbms=csv replace;
run;

