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
  retain transcript_id tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 tpm_bin mean_tpm;
  set tpm_data;
  *if mean_tpm le &tpmQ1. then tpm_bin=1;
  *else if mean_tpm le &tpmQ3. then tpm_bin=2;
  *else tpm_bin=3;
  diff_tpm=(tpm_nsc1-tpm_nsc2);
  mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
  diff_log_tpm=(log_tpm_nsc1-log_tpm_nsc2);
  if mean_log_tpm le 2 then tpm_bin=1;
  else if mean_log_tpm le 4 then tpm_bin=2;
  else tpm_bin=3;
run;


proc sort data=export_data;
   by tpm_bin;
run;


proc corr data=export_data pearson;
  by tpm_bin;
  *var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
  var log_tpm_nsc1 log_tpm_nsc2;
run;

*proc corr data=export_data spearman;
*  by tpm_bin;
*  var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
*run;


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
		N	r2(log)	rho	N	r2	r2(log)	rho	N	r2	r2(log)	rho
PacBio All	119838
RefSeq All	
Event Exp	
Event 100%	
Event 75%	
9Event 75% APN5	
Event 75% APN10	
*/

