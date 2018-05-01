ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Data for example concordance and BA plots for T1D data */

%macro concord(datain,dataout);

data sample1;
   set eventloc.&datain.;
   where library="SL30324";
   keep tpm transcript_id;
   rename tpm=tpm_sample1;
run;

data sample2;
   set eventloc.&datain.;
   where library="SL30327";
   keep tpm transcript_id;
   rename tpm=tpm_sample2;
run;

proc sort data=sample1;
  by transcript_id;
proc sort data=sample2;
  by transcript_id;
run;

data merged;
  merge sample1 (in=in1) sample2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data tpm_data;
  set merged;
  log_tpm_sample1=log(tpm_sample1+1);
  log_tpm_sample2=log(tpm_sample2+1);
  mean_tpm=mean(tpm_sample1,tpm_sample2);
  mean_log_tpm=mean(log_tpm_sample1,log_tpm_sample2);
  diff_log_tpm=log_tpm_sample1-log_tpm_sample2;
run;

proc means data=tpm_data noprint;
  where mean_tpm > 0;
  var mean_tpm;
  output out=tpm_distrib min=min q1=q1 q3=q3 max=max;
run;

%local tpmMin tpmQ1 tpmQ3 tpmMax;

data _null_;
   set tpm_distrib;
   call symputx('tpmMin',min);
   call symputx('tpmQ1',q1);
   call symputx('tpmQ3',q3);
   call symputx('tpmMax',max);
run;

data export_data;
  retain transcript_id tpm_sample1 tpm_sample2 log_tpm_sample1 log_tpm_sample2 tpm_bin mean_tpm;
  set tpm_data;
  if mean_tpm le &tpmQ1. then tpm_bin=1;
  else if mean_tpm le &tpmQ3. then tpm_bin=2;
  else tpm_bin=3;
run;

proc sort data=export_data;
   by tpm_bin;
proc corr data=export_data pearson;
  by tpm_bin;
  var log_tpm_sample1 log_tpm_sample2;
run;



proc sort data=export_data;
   by descending tpm_bin;
run;

proc export data=export_data outfile="!MCLAB/event_analysis/analysis_output/rsem_t1d_concordance_example_&dataout..csv"
   dbms=csv replace;
run;

%mend;

%concord(hg19_rsem_all_xscripts,all);
%concord(hg19_rsem_75perc_apn5_xscripts,reduced);


/* CONCORDANCE

BIN	FULL			REDUCED
	N	R		N	R
Low	378690	0.002	9407	0.053
Med	144426	0.002	11354	0.333
Hi	72971	0.781	5676	0.867

*/

