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

  /* Binning 1: off of lowest TPM */
  if mean_log_tpm=0 then bin1=0;
  else if min_log_tpm < 0.5 then bin1=1;
  else if min_log_tpm < 2.5 then bin1=2;
  else if min_log_tpm < 5 then bin1=3;
  else bin1=4;


  /* Binning 2: off of mean TPM */
  if mean_log_tpm=0 then bin2=0;
    else if mean_log_tpm < 0.5 then bin2=1;
    else if mean_log_tpm < 2.5 then bin2=2;
    else if mean_log_tpm < 5 then bin2=3;
    else bin2=4;

  /* Binning 3: combination */
  if mean_log_tpm=0 then bin3=0;
  else if mean_log_tpm < 0.5 or min_log_tpm < 0.5 then bin3=1;
  else if mean_log_tpm < 2.5 and min_log_tpm < 2.5 then bin3=2;
  else if mean_log_tpm < 5 or min_log_tpm < 5 then bin3=3;
  else if mean_log_tpm ge 5 and min_log_tpm ge 5 then bin3=4;
  else bin3=5;

  /* bin 4: 10 levels */
  if mean_log_tpm=0 then bin4=0;
  else if min_log_tpm < 1 then bin4=1;
  else if min_log_tpm < 3 then bin4=2;
  else if min_log_tpm < 5 then bin4=3;
  else bin4=4;

  diff_tpm=(tpm_nsc1-tpm_nsc2);
  diff_log_tpm=(log_tpm_nsc1-log_tpm_nsc2);
run;


proc sort data=export_data;
   by bin4;
proc corr data=export_data pearson;
  by bin4;
  var tpm_nsc1 tpm_nsc2;
run;   


proc sort data=export_data;
   by descending mean_log_tpm;
run;


data export_data2;
  set export_data;
  if mean_log_tpm=0 then delete;
run;

*proc export data=export_data2 outfile="!MCLAB/event_analysis/analysis_output/rsem_tpm_&datain..csv"
*dbms=csv replace;
*run;

%mend;

%export(pacbio_all);
%export(refseq_all);
%export(events_exp_100perc);
%export(events_exp_75perc_apn5);

data check;
  set export_data2;
  where min_log_tpm < 2;
run;



bin	PB	RS	EA100	EA75
0	.	.	.	.
1	026
2	424
3	
4	938


bin	PB	RS	EA100	EA75
0	.	.	.	.
1	027	002	010	050
2	193	087	208	122
3	940	967	964	964


bin	PB	RS
0	.	.
1	-288	249
2	-147	059
3	-119	001

4	-103	044
5	-026	052
6	065	072
7	095	144
8	144	157
9	142	239
10	949	967


bin	PB	RS	EA100	EA75
0	.	.	.	.
1	027	002	010	059
2	112	000	062	051
3	058	030	076	028
4	630	446	272	211
5	935	964	962	962

bin	PB	RS	EA100	EA75
0	.	.	.	.
1	027	002	010	059
2	112	000	062	051
3	818	837	817	774
4	936	963	961	961

bin	PB	RS	EA100	EA75
0	.	.	.	.
1	027	002	010	059
2	112	000	062	051
3	629	737	698	645
4	939	966	963	963


bin	PB	RS	EA100	EA75
0	.	.	.	.
1	027	002	-010	-059
2	379	572	486	433
3	747	719	724	673
4	936	963	961	961


bin	PB	RS	PB	RS	EA100	EA75
0	.	.	.	.	.	.
1	-091	-002	-027	002	-010	-059	
2	-021	021	-112	000	-062	-050
3	414	170	305	103	124	158
4	939	967	939	967	964	963

proc sort data=export_data;
   by bin1;
proc corr data=export_data pearson;
  by bin1;
  var log_tpm_nsc1 log_tpm_nsc2;
run;   
proc sort data=export_data;
   by bin2;
proc corr data=export_data pearson;
  by bin2;
  var log_tpm_nsc1 log_tpm_nsc2;
run;   

proc sort data=export_data;
   by bin3;
proc corr data=export_data pearson;
  by bin3;
  var log_tpm_nsc1 log_tpm_nsc2;
run;   

  if mean_log_tpm=0 then bin4=0;
  else if min_log_tpm < 0.5 then bin4=1;
  else if min_log_tpm < 1.0 then bin4=2;
  else if min_log_tpm < 1.5 then bin4=3;
  else if min_log_tpm < 2.0 then bin4=4;
  else if min_log_tpm < 2.5 then bin4=5;
  else if min_log_tpm < 3.0 then bin4=6;
  else if min_log_tpm < 3.5 then bin4=7;
  else if min_log_tpm < 4.0 then bin4=8;
  else if min_log_tpm < 4.5 then bin4=9;
  else bin4=10;




data check;
  set export_data;
  where bin3=5;
run;


		Bin0		Bin1		 Bin2		Bin3		Bin 4		
		N	r2(log)	N	r2(log)	 N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	1002	n.d.	1635  	-0.28793 6328	0.71356	6135	0.89128	1004	0.93839	
		1002	n.d.	1070	-0.33702 6273	0.56759	6657	0.85759	1102	0.92904
		1002	n.d.	1635  	-0.28793 5709	0.70760	6754	0.86867	1004	0.93839	
RefSeq (all)	62940	n.d.	45378	0.24970	 14341	0.79972	5009	0.90001	963	0.96401
		62940	n.d.	42903	0.33531	 16275	0.68247	5463	0.83181	1050	0.95180
		62940	n.d.	45378	0.24970	 13803	0.79766	5547	0.85612	963	0.96401	




		Bin0		Bin1		 Bin2		Bin3		Bin 4
		N	r2(log)	N	r2(log)	 N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	1002	n.d.	1635  	-0.28793 6328	0.71356	6135	0.89128	1004	0.93839	
		1002	n.d.	1070	-0.33702 6273	0.56759	6657	0.85759	1102	0.92904
		1002	n.d.	1635	-0.28793 6328	0.71356	6135	0.89128	1004	0.93839
RefSeq (all)	62940	n.d.	45378	0.24970	 14341	0.79972	5009	0.90001	963	0.96401
		62940	n.d.	42903	0.33531	 16275	0.68247	5463	0.83181	1050	0.95180
		62940	n.d.	45378	0.24970	 14341	0.79972	5009	0.90001	963	0.96401
EA 100, apn0	


EA 75, apn5	






/* 
		Bin1		Bin2		Bin3
		N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	6944	0.71871	6649	0.81780	2511	0.95128
RefSeq (all)	120544	0.87317	5940	0.83327	2147	0.96964
EA 100, apn0	8773	0.79400	4032	0.81770	1929	0.96917
EA 75, apn5	7878	0.77202	3925	0.80748	1937	0.97048
0.34199
		Bin0		Bin1		Bin2		Bin3		Bin 4
		N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	1002	n.d.	2727	-0.0057	3215	0.28834	6649	0.81780	2511	0.95128
RefSeq (all)	62940	n.d.	51197	0.60198	6407	0.42881	5940	0.83327	2147	0.96964
EA 100, apn0	970	n.d.	4975	0.34199	2828	0.37750	4032	0.81770	1929	0.96917
EA 75, apn5	1868	n.d.	3964	0.20251	2086	0.28536	3925	0.80748	1937	0.97048


		Bin1		Bin2		Bin3		Bin 4
		N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	3729	0.26541	3215	0.28834	6649	0.81780	2511	0.95128
RefSeq (all)	114137	0.68206	6407	0.42881	5940	0.83327	2147	0.96964
EA 100, apn0	5945	0.44785	2828	0.37750	4032	0.81770	1929	0.96917
EA 75, apn5	5792	0.41789	2086	0.28536	3925	0.80748	1937	0.97048


12226
		Bin0		Bin1		Bin2		Bin3		Bin 4
		N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	1002	n.d.	1635   -0.02727	4307	0.19309	6649	0.63038	2511	0.93588
RefSeq (all)	62940	n.d.	45378	0.00290	12226	0.08763	5940	0.44637	2147	0.96455
EA 100, apn0	970	n.d.	3377   -0.01020	4426	0.20866	4032	0.27275	1929	0.96203
EA 75, apn5	1868	n.d.	2715   -0.05926	3295	0.12282	3925	0.21152	1973	0.96259


		N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	N	r2(log)	
PacBio (all)	1002	n.d.	1635   -0.28793	4307	0.53424	6649	0.81780	2511	0.95128
RefSeq (all)	62940	n.d.	45378	0.24970	12226	0.68269	5940	0.83327	2147	0.96964
EA 100, apn0	970	n.d.	3377   -0.05063	4426	0.62850	4032	0.81770	1929	0.96917
EA 75, apn5	1868	n.d.	2715   -0.16567	3295	0.56251	3925	0.80748	1937	0.97048


Try:
<0.5,

*/

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
