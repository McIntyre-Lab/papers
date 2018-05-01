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
   keep transcript_id  mean_tpm tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
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
  set tpm_data;
  if mean_tpm le &tpmQ1. then tpm_bin=1;
  else if mean_tpm le &tpmQ3. then tpm_bin=2;
  else tpm_bin=3;
  drop mean_tpm;
run;

proc sort data=export_data;
   by tpm_bin;
run;


proc corr data=export_data pearson;
  by tpm_bin;
  var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
run;

proc corr data=export_data spearman;
  by tpm_bin;
  var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
run;


   proc export data=export_data outfile="!MCLAB/event_analysis/analysis_output/rsem_fpkm_&datain..csv"
   dbms=csv replace;
run;

%mend;


%export(refseq_all);
%export(pacbio_all);
%export(events_exp_any);
%export(events_exp_100perc);
%export(events_exp_75perc);
%export(events_exp_75perc_apn5);
%export(events_exp_75perc_apn10);


/* CORR BY BIN
		Bin1				Bin2				Bin3
		N	r2	r2(log)	rho	N	r2	r2(log)	rho	N	r2	r2(log)	rho
RefSeq All	80455	0.11126	0.11427	0.28539	31757	0.49074	0.53274	0.40681	16419	0.96805	0.94075	0.89585	
Event Exp	39553	0.07968	0.08290	0.25646	22618	0.52411	0.57104	0.44153	11364	0.96581	0.94602	0.90768
Event 75%	14710	0.15433	0.18013	0.37805	13258	0.66074	0.66462	0.65236	6654	0.96489	0.94640	0.92001
Event 75% APN5	4841	0.38807	0.48099	0.58100	5931	0.79641	0.76827	0.82718	2968	0.96243	0.95035	0.92521
Event 100%	4414	0.39839	0.48106	0.57050	6880	0.79014	0.78045	0.82303	3440	0.96296	0.94892	0.92393
Event 75% APN10	3335	0.44808	0.56321	0.63960	4363	0.83023	0.79988	0.86259	2181	0.96557	0.94797	0.92108
PacBio All	4779	0.53537	0.66089	0.68146	7550	0.73962	0.68778	0.78093	3775	0.93748	0.92512	0.88866
*/

