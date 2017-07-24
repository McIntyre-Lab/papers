ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";


/* BA plots for each RSEM output. I am going to calculate the data needed for plots here then export
   to python for plotting

   Put in macro so I can iterate through the data

*/

%macro baPlot(datain);

data calc_mean_diff;
   set event.rsem_&datain.;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
   diff_tpm=(tpm_nsc1-tpm_nsc2);

   mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
   diff_log_tpm=(log_tpm_nsc1-log_tpm_nsc2);

   keep transcript_id mean_tpm diff_tpm mean_log_tpm diff_log_tpm;
run;

   proc export data=calc_mean_diff outfile="!MCLAB/event_analysis/analysis_output/rsem_ba_plot_data_&datain..csv"
   dbms=csv replace;
run;

%mend;


%baPlot(refseq_all);
%baPlot(pacbio_all);
%baPlot(events_exp_any);
%baPlot(events_exp_100perc);
%baPlot(events_exp_75perc);
%baPlot(events_exp_75perc_apn5);
%baPlot(events_exp_75perc_apn10);


