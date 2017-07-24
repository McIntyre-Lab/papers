ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/*
Want to plot the mean TPM against the % total events detected
Going to use the "all Event Analysis transcripts" TPM for this.
*/

data xs_perc_dtct;
   set event.xscripts_w_unique_by_bin;
   keep transcript_id perc_features_dtct;
run;

data xs_tpm;
   set event.rsem_events_exp_any;
   mean_tpm=(log(tpm_nsc1+1)+log(tpm_nsc2+1))/2;
   keep transcript_id mean_tpm;
run;

proc sort data=xs_perc_dtct;
   by transcript_id;
proc sort data=xs_tpm;
   by transcript_id;
run;

data xs_tpm2dtct;
  merge xs_perc_dtct (in=in1) xs_tpm (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* sort by X-values */

proc sort data=xs_tpm2dtct;
   by perc_features_dtct;
run;


proc export data=xs_tpm2dtct
     outfile="!MCLAB/event_analysis/analysis_output/xscripts_mean_tpm2perc_dtct.csv"
     dbms=csv replace;
run;

