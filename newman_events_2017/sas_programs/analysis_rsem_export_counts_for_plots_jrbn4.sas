ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";


/* Export RSEM FPKM data for making plots -- I also want to bin transcripts into the same 3 bins as
  in the SQANTI paper: min-Q1, Q1-Q3, Q3-max */

%macro export(datain);

data tpm_data;
   set event.rsem_&datain.;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_tpm tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
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
  retain transcript_id tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 tpm_bin mean_tpm;
  set tpm_data;
  if mean_tpm le &tpmQ1. then tpm_bin=1;
  else if mean_tpm le &tpmQ3. then tpm_bin=2;
  else tpm_bin=3;
  diff_tpm=(tpm_nsc1-tpm_nsc2);
  mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
  diff_log_tpm=(log_tpm_nsc1-log_tpm_nsc2);
run;


proc sort data=export_data;
   by descending tpm_bin;
*proc corr data=export_data pearson;
*  by tpm_bin;
*  *var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
*  var log_tpm_nsc1 log_tpm_nsc2;
run;

   proc export data=export_data outfile="!MCLAB/event_analysis/analysis_output/rsem_fpkm_&datain..csv"
   dbms=csv replace;
run;

%mend;

%export(pacbio_all);
%export(refseq_all);
%export(events_exp_any);
%export(events_exp_100perc);
%export(events_exp_75perc);
%export(events_exp_75perc_apn5);
%export(events_exp_75perc_apn10);
%export(events_exp_100perc_apn5);
%export(events_exp_100perc_apn10);
%export(events_exp_50perc);
%export(events_exp_50perc_apn5);
%export(events_exp_50perc_apn10);
%export(pacbio_apn0);
%export(pacbio_apn5);
%export(pacbio_apn10);


/* CORR BY BIN:  TPM low-med-high

		Bin1		Bin2		Bin3
		N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	4779	0.66089	7550	0.68778	3775	0.92512
RefSeq (all)	80455	0.11427	31757	0.53274	16419	0.94075
Events (any)	39553	0.08290	22618	0.57104	11364	0.94602
EA 100, apn0	4414	0.48106	6880	0.78045	3440	0.94892
EA 75, apn0	14710	0.18013	13258	0.66462	6654	0.94640
EA 75, apn5	4841	0.48099	5931	0.76827	2968	0.95035
EA 75, apn10	3335	0.56321	4363	0.79988	2181	0.94797
EA 100, apn5	1031	0.84396	1856	0.82083	928	0.92661
EA 100, apn10	698	0.89560	1272	0.81554	636	0.90741
EA 50, apn0	20955	0.16542	16553	0.63315	8375	0.94530
EA 50, apn5	8000	0.32526	8219	0.73486	4117	0.94594
EA 50, apn10	5567	0.42273	6356	0.78321	3179	0.95234
PacBio, apn0	1703	0.83907	3056	0.80326	1527	0.94130
PacBio, apn5	1646	0.83298	2955	0.79660	1477	0.93997
PacBio, apn10	1599	0.84224 2857	0.79688	1428	0.94189
*/
