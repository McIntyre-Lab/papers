ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";


/* Export RSEM FPKM data for making plots -- I also want to bin transcripts into low/low-med/med-high/high exp bins:
   bin 1 = lowest is <1
   bin 2 = lowest is between 1 and 2
   bin 3 = lowest is between 2 and 4
   bin 4 = lowest is >4
   */

%macro kappa(datain);

data tpm_data;
   set event.rsem_&datain.;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
   min_log_tpm=min(log_tpm_nsc1,log_tpm_nsc2);
   mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
   min_tpm=min(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 min_tpm min_log_tpm mean_log_tpm ;
run;


data agreement_data;
  set tpm_data;
  /* bin 4 levels */
  if mean_log_tpm=0 then tpm_bin=0;
  else if min_log_tpm < 0.5 then tpm_bin=1;
  else if min_log_tpm < 2 then tpm_bin=2;
  else if min_log_tpm < 4 then tpm_bin=3;
  else tpm_bin=4;

  /* Set agreements */
  if tpm_nsc1=0 then flag_nsc1_tpm_gt0=0;
  else flag_nsc1_tpm_gt0=1;

  if tpm_nsc2=0 then flag_nsc2_tpm_gt0=0;
  else flag_nsc2_tpm_gt0=1;
run;

data agreement_data2;
   set agreement_Data;
   if flag_nsc1_tpm_gt0=0 and flag_nsc2_tpm_gt0=0 then delete;
run;


proc freq data=agreement_data2 noprint;
   where flag_nsc1_tpm_gt0 ne 0 and flag_nsc2_tpm_gt0 ne 0;
   tables flag_nsc1_tpm_gt0*flag_nsc2_tpm_gt0 / out=counts;
   test agree;
   output out=kappa_stats AGREE;
run;

data counts2;
   length category $8.;
   set counts;
   if flag_nsc1_tpm_gt0=0 and flag_nsc2_tpm_gt0=0 then category="N_0v0";
   else if flag_nsc1_tpm_gt0=0 and flag_nsc2_tpm_gt0=1 then category="N_0v1";
   else if flag_nsc1_tpm_gt0=1 and flag_nsc2_tpm_gt0=0 then category="N_1v0";
   else if flag_nsc1_tpm_gt0=1 and flag_nsc2_tpm_gt0=1 then category="N_1v1";
   keep category count;
run;

proc transpose data=counts2 out=counts_sbys;
   id category;
   var count;
run;

data freq_&datain;
   length xscript_set $30.;
   set counts_sbys;
   xscript_set="&datain.";
   perc_disagree=(N_0v1+N_1v0)/(N_0v1+N_1v0+N_0v0+N_1v1)*100;
   drop _NAME_ _LABEL_;
run;
  

data agree_&datain.;
   length xscript_set $30.;
   set kappa_stats;
   xscript_set="&datain.";
run;

%mend;

%kappa(pacbio_all);
%kappa(refseq_all);
%kappa(events_exp_100perc);
%kappa(events_exp_75perc_apn5);
%kappa(events_exp_any);
%kappa(events_exp_75perc);
%kappa(events_exp_100perc_apn5);
%kappa(events_exp_75perc_apn10);
%kappa(events_exp_100perc_apn10);
%kappa(events_exp_50perc);
%kappa(events_exp_50perc_apn5);
%kappa(events_exp_50perc_apn10);
%kappa(pacbio_apn0);
%kappa(pacbio_apn5);
%kappa(pacbio_apn10);

data all_kappa_data;
  set agree_: ;
run;


data all_freq_data;
  set freq_: ;
run;

proc sort data=all_kappa_data;
   by xscript_set;
proc sort data=all_freq_data;
   by xscript_set;
run;

data all_agreement_data;
  merge all_freq_data (in=in1) all_kappa_data (in=in2);
  by xscript_set;
  if in1 and in2;
run;

proc export data=all_agreement_data
     outfile="!MCLAB/event_analysis/analysis_output/agreement_stats_by_transcript_set.csv"
     dbms=csv replace;
run;


/*

                                                                 perc_
Obs   xscript_set              N_0v0   N_0v1   N_1v0   N_1v1   disagree      N     _MCNEM_   P_MCNEM

 1    events_exp_100perc         970     532     563   12669     7.4318    14734    0.8776   0.34885
 2    events_exp_75perc_apn5    1868     439     757   10676     8.7045    13740   84.5518   0.00000
 3    pacbio_all                1002     303     508   14291     5.0360    16104   51.8187   0.00000
 4    refseq_all               62940   11528   10441   43722    17.0791   128631   53.7835   0.00000

Obs   _KAPPA_    E_KAPPA   L_KAPPA   U_KAPPA     E0_KAPPA   Z_KAPPA   PL_KAPPA   PR_KAPPA   P2_KAPPA

 1    0.59779   0.011031   0.57617   0.61941   .008237805    72.566       .          0          0
 2    0.70473   0.007982   0.68908   0.72037   .008504799    82.862       .          0          0
 3    0.68447   0.010395   0.66410   0.70484   .007855021    87.138       .          0          0
 4    0.65065   0.002142   0.64645   0.65484   .002787804   233.390       .          0          0

*/

