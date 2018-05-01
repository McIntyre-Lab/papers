ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';


/* For each set of transcripts, look at simple agreement stats for eXpress data, as well as the by-tpm_bin correlation */

%macro kappa(outName);

data NSC1;
   set event.xprs_NSC1_&outName.;
   keep target_id tpm;
   rename target_id=transcript_id tpm=tpm_nsc1;
run;

data NSC2;
   set event.xprs_NSC2_&outName.;
   keep target_id tpm;
   rename target_id=transcript_id tpm=tpm_nsc2;
run;

proc sort data=NSC1;
  by transcript_id;
proc sort data=NSC2;
  by transcript_id;
run;

data tpm_data;
  merge NSC1 (in=in1) NSC2 (in=in2);
  by transcript_id;
  if in1 and in2;
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
   *if flag_nsc1_tpm_gt0=0 and flag_nsc2_tpm_gt0=0 then delete;
run;

/* Calc by-bin correlation */

proc sort data=agreement_data;
  by tpm_bin;
proc corr data=agreement_data pearson;
  by tpm_bin;
  var log_tpm_nsc1 log_tpm_nsc2;
run;

proc freq data=agreement_data2 noprint;
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

data freq_&outName;
   retain xscript_set perc_disagree N_0v0 N_0v1 N_1v0 N_1v1;
   format N_0v0 best12.;
   format N_0v1 best12.;
   format N_1v0 best12.;
   format N_1v1 best12.;
   length xscript_set $30.;
   set counts_sbys;
   xscript_set="&outName.";
   if N_0v0=. then N_0v0=0;
   if N_0v1=. then N_0v1=0;
   if N_1v0=. then N_1v0=0;
   if N_1v1=. then N_1v1=0;
   perc_disagree=(N_0v1+N_1v0)/(N_0v1+N_1v0+N_0v0+N_1v1)*100;
   drop _NAME_ _LABEL_;
run;
  
  
data agree_&outName.;
   length xscript_set $30.;
   set kappa_stats;
   xscript_set="&outName.";
run;

%mend;

%kappa(refseq_all);
%kappa(events_100prc_apn0);
%kappa(events_75prc_apn5);

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
  if sum(N_0v0,N_0v1,N_1v0, N_1v1) ne N then do;
      *if these don't agree then kappa stats weren't calculated so set all to missing;
         _NCNEM_=.; P_MCNEM=.; _KAPPA_=.; E_KAPPA=.; L_KAPPA=.; U_KAPPA=.; E0_KAPPA=.; Z_KAPPA=.;
    end;
run;

proc export data=all_agreement_data
     outfile="!MCLAB/event_analysis/analysis_output/agreement_stats_by_transcript_set_express.csv"
     dbms=csv replace;
run;

/*
		bin0		bin1		 bin2		bin3		bin4
		N	r	N	r	 N	r	N	r	N	r	
Refseq All	31289	n.d	64319	0.78725	 22314	0.86802	8349	0.89852	2360	0.97283
EA 100% APN0	0	n.d	146	-0.34217 5553	0.82613	5974	0.88409	2067	0.97509
EA 75% APN5	117	n.d	1735	0.36625	 5128	0.79127	5648	0.86996	2106	0.97487


*/

