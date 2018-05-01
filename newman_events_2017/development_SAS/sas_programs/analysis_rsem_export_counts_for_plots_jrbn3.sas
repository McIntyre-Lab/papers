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

data xs2dtct;
   set event.xscripts_w_unique_by_bin_total;
   keep transcript_id perc_features_dtct;
run;

data xs2type;
   set event.xscripts_w_unique_by_bin;
   keep transcript_id flag_xscript_has_unique flag_xscript_has_unique_dtct;
run;

proc sort data=xs2dtct;
   by transcript_id;
proc sort data=xs2type;
   by transcript_id;
run;



%if &datain.=pacbio_all or &datain.=pacbio_apn0 or &datain.=pacbio_apn5 or &datain.=pacbio_apn10
   %then %do;
data pb2xs;
  set event.pacbio2refseq_xscript_nomulti;
  keep pacbio_id transcript_id;
run;

proc sort data=pb2xs;
   by transcript_id pacbio_id;
run;

data pb2xs_w_type;
   merge pb2xs (in=in1) xs2dtct (in=in2) xs2type (in=in3);
   by transcript_id;
   if not in2 then perc_features_dtct=-9;
   if not in3 then do; flag_xscript_has_unique=-9; flag_xscript_has_unique_dtct=-9; end;
   if in1 then output;
   drop transcript_id;
   rename pacbio_id=transcript_id;
run;

proc sort data=pb2xs_w_type;
   by transcript_id;
proc sort data=tpm_data;
   by transcript_id;
run;

data tpm_data2;
   merge tpm_data (in=in1) pb2xs_w_type (in=in2);
   by transcript_id;
   if in1 then output;
run;

%end;
%else %do;

proc sort data=tpm_data;
   by transcript_id;
run;

data tpm_data2;
   merge tpm_data (in=in1) xs2dtct (in=in2) xs2type (in=in3);
   by transcript_id;
   if not in2 then perc_features_dtct=-9;
   if not in3 then do; flag_xscript_has_unique=-9; flag_xscript_has_unique_dtct=-9; end;
   if in1 then output;
run;

%end;

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
  set tpm_data2;
  if mean_tpm le &tpmQ1. then tpm_bin=1;
  else if mean_tpm le &tpmQ3. then tpm_bin=2;
  else tpm_bin=3;
  diff_tpm=(tpm_nsc1-tpm_nsc2);
  mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
  diff_log_tpm=(log_tpm_nsc1-log_tpm_nsc2);
  /* detection bin */
  if perc_features_dtct=1 then dtct_bin=1; *100% features on;
  else if  perc_features_dtct ge 0.75 then dtct_bin=2; *75% features on;
  else if  perc_features_dtct ge 0.5 then dtct_bin=3; *50% features on;
  else if  perc_features_dtct ge 0.25 then dtct_bin=4; *25% features on;
  else if  perc_features_dtct > 0 then dtct_bin=5; *1% features on;
  else if  perc_features_dtct = 0 then dtct_bin=6; *All features off;
  else dtct_bin=7; *not included in analysis;
   /* uniqueness bin */
  if flag_xscript_has_unique=1 and flag_xscript_has_unique_dtct=1 then uniq_bin=1; *has uniq detected;
  else if flag_xscript_has_unique=0 and flag_xscript_has_unique_dtct=1 then uniq_bin=-9; *should be empty;
  else if flag_xscript_has_unique=1 and flag_xscript_has_unique_dtct=0 then uniq_bin=2; *has uniq no dtct;
  else if flag_xscript_has_unique=0 and flag_xscript_has_unique_dtct=0 then uniq_bin=3; *no uniq;
  else uniq_bin=4; *excluded;
run;


proc sort data=export_data;
   by tpm_bin;
proc corr data=export_data pearson;
  by tpm_bin;
  *var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
  var log_tpm_nsc1 log_tpm_nsc2;
run;

proc sort data=export_data;
   by dtct_bin;
proc corr data=export_data pearson;
  by dtct_bin;
  *var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
  var log_tpm_nsc1 log_tpm_nsc2;
run;

proc sort data=export_data;
   by uniq_bin;
proc corr data=export_data pearson;
  by uniq_bin;
  *var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
  var log_tpm_nsc1 log_tpm_nsc2;
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


   CORR BY BIN: Perc dtct

		Bin1		Bin2		Bin3		Bin4		Bin5		Bin6		Bin7
		N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)
PacBio (all)	4877	0.96362	1215	0.94818	144	0.97666	26	0.97087	24	0.88548	0	n.d	9818	0.95187
RefSeq (all)	14734	0.97172	19888	0.93945	11261	0.92150	11689	0.90411	15963	0.92133	0	n.d	55096	0.97069
Events (any)	14734	0.97114	19888	0.93744	11261	0.93394	11689	0.94977	15963	0.96027
EA 100, apn0	14734	0.97238
EA 75, apn0	14734	0.97064	19888	0.93739
EA 75, apn5	9310	0.96990	4430	0.95327	
EA 75, apn10	7259	0.97182	2620	0.95818
EA 100, apn5	3815	0.96942	
EA 100, apn10	2606	0.97094
EA 50, apn0	14734	0.97080	19888	0.93702	11261	0.93490
EA 50, apn5	11256	0.96878	8805	0.94451	275	0.95573
EA 50, apn10	9572	0.97132	5363	0.95282	167	0.97349
PacBio, apn0	4877	0.96701	1215	0.95817	144	0.97234	26	0.97892	24	0.88623
PacBio, apn5	4850	0.96691	1159	0.95692	68	0.95167	1	n.d		
PacBio, apn10	4775	0.96671	1064	0.95863	44	0.96407	1	n.d


   CORR BY BIN: uniqueness

		Bin1		Bin2		Bin3		Bin4	
		N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	4084	0.96936	312	0.96649	1890	0.95067	9818	0.95187
RefSeq (all)	22319	0.97812	16316	0.92790	34900	0.95797	55096	0.97069
Events (any)	22319	0.97736	16316	0.92736	34900	0.96106
EA 100, apn0	9730	0.97505	0	n.d	5004	0.96046
EA 75, apn0	16412	0.97414	5117	0.92195	13093	0.95404
EA 75, apn5	7958	0.97577	1230	0.93442	4552	0.95506
EA 75, apn10	6003	0.97774	733	0.93471	3143	0.95745
EA 100, apn5	2762	0.96849	0	n.d	1053	0.96357
EA 100, apn10	1950	0.96936	0	n.d	656	0.96655
EA 50, apn0	18838	0.97531	8115	0.92083	18930	0.95645
EA 50, apn5	10713	0.97437	2527	0.92544	7096	0.95360
EA 50, apn10	8578	0.97709	1553	0.93258	4971	0.95710
PacBio, apn0	4084	0.97068	312	0.97142	1890	0.96248
PacBio, apn5	3982	0.96881	263	0.96308	1833	0.95994
PacBio, apn10	3858	0.96857	237	0.96490	1789	0.95955
*/


