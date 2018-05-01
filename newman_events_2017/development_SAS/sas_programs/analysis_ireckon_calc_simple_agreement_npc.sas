ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';


/* For each set of transcripts, look at simple agreement stats for eXpress data, as well as the by-tpm_bin correlation */

data NSC1 NSC2;
   set event.ireckon_results_nsc;
   if sample_id="NSC1" then output NSC1;
   if sample_id="NSC2" then output NSC2;
   drop rpkm_bin;
run;

data NSC1_2;
  set NSC1;
  rename rpkm=rpkm_nsc1 log_rpkm=log_rpkm_nsc1;
run;

data NSC2_2;
  set NSC2;
  rename rpkm=rpkm_nsc2 log_rpkm=log_rpkm_nsc2;
run;

/* where there are multiple entries for the same transcript, take the larger value */


proc sort data=NSC1_2;
  by chr start stop strand gene_id transcript_id ;
proc sort data=NSC2_2;
  by chr start stop strand  gene_id transcript_id;
run;

proc means data=NSC1_2 noprint;
  by chr start stop strand gene_id transcript_id;
  var rpkm_nsc1 log_rpkm_nsc1;
  output out=NSC1_3(drop=_TYPE_ _FREQ_) max=;
run;

proc means data=NSC2_2 noprint;
  by chr start stop strand gene_id transcript_id;
  var rpkm_nsc2 log_rpkm_nsc2;
  output out=NSC2_3(drop=_TYPE_ _FREQ_) max=;
run;


data rpkm_data;
  merge NSC1_3 (in=in1) NSC2_3 (in=in2);
  by chr start stop strand  gene_id transcript_id;
  if rpkm_nsc1=. then rpkm_nsc1=0;
  if rpkm_nsc2=. then rpkm_nsc2=0;
  if log_rpkm_nsc1=. then log_rpkm_nsc1=0;
  if log_rpkm_nsc2=. then log_rpkm_nsc2=0;
  mean_rpkm=(rpkm_nsc1+rpkm_nsc2)/2;
  min_log_rpkm=min(log_rpkm_nsc1,log_rpkm_nsc2);
  mean_log_rpkm=(log_rpkm_nsc1+log_rpkm_nsc2)/2;
  min_rpkm=min(rpkm_nsc1,rpkm_nsc2);
  keep chr start stop strand gene_id transcript_id mean_rpkm rpkm_nsc1
       rpkm_nsc2 log_rpkm_nsc1 log_rpkm_nsc2
       min_rpkm min_log_rpkm mean_log_rpkm ;
run;


proc univariate data=rpkm_data;
    where min_log_rpkm > 0;
    var min_log_rpkm;
run;

/*


Quantile          Estimate

100% Max       8.37088E+00
99%            4.25692E+00
95%            2.45097E+00
90%            1.52554E+00
75% Q3         6.31921E-01
50% Median     2.05537E-01
25% Q1         6.88595E-02
10%            1.86357E-02
5%             5.48884E-04
1%             3.69653E-05
0% Min         1.04803E-05

My thresholds for RSEM/eXpress for bins are:
Bin	logTPM	Approximate range of distribution
0 	0	Minimum
1 	0-0.5	Min -> Q1
2	0.5-2	Q1 -> Q3
3	2-4	Q3->90%
4	4-Max	90%->100%

Therefore, for iReckon, my logRPKM bins are

Bin	logRPKM	Range
0	0	Minimum
1	0-0.1	Min -> Q1
2	0.1-0.5	Q1 -> Q3
3	0.5-1.5	Q3 -> 90%
4	1.5+	90% -> 100%

*/





data agreement_data;
  set rpkm_data;
  /* bin 4 levels */
  if mean_log_rpkm=0 then rpkm_bin=0;
  else if min_log_rpkm < 0.1 then rpkm_bin=1;
  else if min_log_rpkm < 0.5 then rpkm_bin=2;
  else if min_log_rpkm < 1.5 then rpkm_bin=3;
  else rpkm_bin=4;

  /* Set agreements */
  if rpkm_nsc1=0 then flag_nsc1_rpkm_gt0=0;
  else flag_nsc1_rpkm_gt0=1;

  if rpkm_nsc2=0 then flag_nsc2_rpkm_gt0=0;
  else flag_nsc2_rpkm_gt0=1;
run;


data agreement_data2;
   set agreement_Data;
   *if flag_nsc1_rpkm_gt0=0 and flag_nsc2_rpkm_gt0=0 then delete;
run;

/* Calc by-bin correlation */

proc sort data=agreement_data;
  by rpkm_bin;
proc corr data=agreement_data pearson;
  by rpkm_bin;
  var log_rpkm_nsc1 log_rpkm_nsc2;
run;

proc freq data=agreement_data2 noprint;
   tables flag_nsc1_rpkm_gt0*flag_nsc2_rpkm_gt0 / out=counts;
   test agree;
   output out=kappa_stats AGREE;
run;

data counts2;
   length category $8.;
   set counts;
   if flag_nsc1_rpkm_gt0=0 and flag_nsc2_rpkm_gt0=0 then category="N_0v0";
   else if flag_nsc1_rpkm_gt0=0 and flag_nsc2_rpkm_gt0=1 then category="N_0v1";
   else if flag_nsc1_rpkm_gt0=1 and flag_nsc2_rpkm_gt0=0 then category="N_1v0";
   else if flag_nsc1_rpkm_gt0=1 and flag_nsc2_rpkm_gt0=1 then category="N_1v1";
   keep category count;
run;


proc transpose data=counts2 out=counts_sbys;
   id category;
   var count;
run;

data freq_iReckon;
   retain xscript_set perc_disagree N_0v0 N_0v1 N_1v0 N_1v1;
   format N_0v0 best12.;
   format N_0v1 best12.;
   format N_1v0 best12.;
   format N_1v1 best12.;
   length xscript_set $30.;
   set counts_sbys;
   xscript_set="iReckon";
   if N_0v0=. then N_0v0=0;
   if N_0v1=. then N_0v1=0;
   if N_1v0=. then N_1v0=0;
   if N_1v1=. then N_1v1=0;
   perc_disagree=(N_0v1+N_1v0)/(N_0v1+N_1v0+N_0v0+N_1v1)*100;
   drop _NAME_ _LABEL_;
run;
  
  
data agree_ireckon;
   length xscript_set $30.;
   set kappa_stats;
   xscript_set="iReckon";
run;


proc sort data=agree_ireckon;
   by xscript_set;
proc sort data=freq_iReckon;
   by xscript_set;
run;

data all_agreement_data;
  merge freq_iReckon (in=in1) agree_ireckon (in=in2);
  by xscript_set;
  if in1 and in2;
run;

proc export data=all_agreement_data
     outfile="!MCLAB/event_analysis/analysis_output/agreement_stats_iReckon.csv"
     dbms=csv replace;
run;

/*
		bin0		bin1		 bin2		bin3		bin4
		N	r	N	r	 N	r	N	r	N	r	
iReckon		0	n.d.	22318	-0.066	 859	0.0716	440	0.31837	240	0.89454

*/

