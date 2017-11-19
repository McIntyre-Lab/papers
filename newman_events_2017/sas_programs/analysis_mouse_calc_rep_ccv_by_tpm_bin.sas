ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";


/* Export RSEM FPKM data for making plots -- I also want to bin transcripts into low/low-med/med-high/high exp bins:
   bin 1 = lowest is <1
   bin 2 = lowest is between 1 and 2
   bin 3 = lowest is between 2 and 4
   bin 4 = lowest is >4
   */

%macro export(datain);

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


data export_data;
  set tpm_data;

  /* bin 4 levels */
  if mean_log_tpm=0 then tpm_bin=0;
  else if min_log_tpm < 0.5 then tpm_bin=1;
  else if min_log_tpm < 2 then tpm_bin=2;
  else if min_log_tpm < 4 then tpm_bin=3;
  else tpm_bin=4;

  diff_tpm=abs(tpm_nsc1-tpm_nsc2);
  diff_log_tpm=abs(log_tpm_nsc1-log_tpm_nsc2);
run;

proc sort data=export_data;
  by tpm_bin;
proc means data=export_data noprint;
  by tpm_bin;
  var mean_log_tpm;
  output out=var_stats_by_bin
         var(mean_log_tpm)=var_mean_log  cv(mean_log_tpm)=cv_mean_log;
run;

data var_stats2;
   length bin_N $15.;
   length bin_CV $15.;
   set var_Stats_by_bin;
   bin_N = catt("bin_",tpm_bin,"_N");
   bin_CV = catt("bin_",tpm_bin,"_CV");
   keep bin_N bin_CV _FREQ_ CV_mean_log ;
run;

proc transpose data=var_stats2 out=counts_sbys;
   id bin_N;
   var _FREQ_;
run;

proc transpose data=var_stats2 out=CV_sbys;
   id bin_CV;
   var CV_mean_log;
run;


data counts_sbys2;
   length xscript_set $12.;
   set counts_sbys;
   xscript_set="&datain.";
   keep xscript_set bin_0_N bin_1_N bin_2_N bin_3_N bin_4_N;
run;

data CV_sbys2;
   length xscript_set $12.;
   set CV_sbys;
   xscript_set="&datain.";
   keep xscript_set bin_0_CV bin_1_CV bin_2_CV bin_3_CV bin_4_CV;
run;

data binstats_&datain.;
  merge counts_sbys2 CV_sbys2;
run;

%mend;

%export(pacbio_all);
%export(refseq_all);
%export(events_exp_100perc);
%export(events_exp_75perc_apn5);
%export(events_exp_any);
%export(events_exp_75perc);
%export(events_exp_100perc_apn5);
%export(events_exp_75perc_apn10);
%export(events_exp_100perc_apn10);
%export(events_exp_50perc);
%export(events_exp_50perc_apn5);
%export(events_exp_50perc_apn10);
%export(pacbio_apn0);
%export(pacbio_apn5);
%export(pacbio_apn10);


data all_binstats;
  set binstats_: ;
run;



proc export data=all_binstats
     outfile="!MCLAB/event_analysis/analysis_output/CV_by_TPM_bin_by_transcript_set.csv"
     dbms=csv replace;
run;


/*
	PACBIO		REFSEQ_ALL		EA_100_APN0	EA_75_APN5
BIN	N	CV_MEAN	N	CV_MEAN	N	CV_MEAN	N	CV_MEAN
0	1002	n/a	62940	n/a	970	n/a	1868	n/a
1	1635	76.3796	45378	128.767	3377	87.1348	2715	81.0323
2	4307	29.3791	12226	36.347 	4426	33.0444	3295	33.1822
3	6649	18.4165	5940	19.045	4032	18.5452	3925	18.2734
4	2511	17.7129	2147	18.859	1929	19.7183	1937	19.8243


Bins:
0 = average TPM is 0 (no expression)
1 = minimum TPM is <0.5	(very low expression)
2 = minimum TPM is <2	(low expression)
3 = minimum TPM is <4	(moderate expression)
4 = minimum TPM is ge 4	(high expression)
*/



