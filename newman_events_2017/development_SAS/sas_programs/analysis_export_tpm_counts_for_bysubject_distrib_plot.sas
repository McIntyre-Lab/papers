ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Export log TPM for plotting. I want to check that the distribution of TPM estimates between subjects is similar,
   otherwise may need to normalize

   Need: Name transcript_id cell_type tpm log_tpm
*/

data subject_list;
   set eventloc.subjects_for_rsem_test;
run;

data design;
   set con.design_by_subject_new;
   keep subject_id cell_type name library;
run;

proc sort data=subject_list;
   by subject_id;
proc sort data=design;
   by subject_id;
run;

data design2;
   merge design (in=in1) subject_list (in=in2);
   by subject_id;
   if in1 and in2;
run;

data counts_for_plot_all;
  set eventloc.hg19_rsem_all_xscripts;
  log_tpm=log(tpm+1);
run;

data counts_for_plot_filt;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
  log_tpm=log(tpm+1);
run;

proc sort data=design2;
   by library;
proc sort data=counts_for_plot_all;
   by library;
proc sort data=counts_for_plot_filt;
   by library;
run;

data counts_for_plot_all2;
   merge design2 (in=in1) counts_for_plot_all (in=in2);
   by library;
   if in1 and in2;
   if log_tpm=0 then delete;
run;

data counts_for_plot_filt2;
   merge design2 (in=in1) counts_for_plot_filt (in=in2);
   by library;
   if in1 and in2;
   if log_tpm=0 then delete;
run;

proc export data=counts_for_plot_all2 outfile="/mnt/store/event_sandbox/tpm_by_subject_cell_all.csv"
  dbms=csv replace;
run;

proc export data=counts_for_plot_filt2 outfile="/mnt/store/event_sandbox/tpm_by_subject_cell_filtered.csv"
  dbms=csv replace;
run;

