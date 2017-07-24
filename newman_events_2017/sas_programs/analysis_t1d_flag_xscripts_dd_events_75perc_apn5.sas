ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Trying something else:
   (1) Flag transcripts that are DD
	- (1) not in all cell types for 75% detection at APN>5
	- (2) TPM >= 1, >= 5
   (2) Calculate variances on:
	(1) Transcripts not DD (event detection, TPM1, TPM5)
	(2) Gene-level (all transcripts)
	(3) Gene-level (non-DD transcripts)

*/

/* Flag differentially-detected transcripts 
   If transcript is on in less than all 3 cell types then it is DD */

data flag_xs_dd_event;
   set event.hg19_flag_xscripts_w_unique;
   if perc_features_dtct_cd4 >= 0.75 then flag_cd4_75perc_apn5=1; else flag_cd4_75perc_apn5=0;
   if perc_features_dtct_cd8 >= 0.75 then flag_cd8_75perc_apn5=1; else flag_cd8_75perc_apn5=0;
   if perc_features_dtct_cd19 >= 0.75 then flag_cd19_75perc_apn5=1; else flag_cd19_75perc_apn5=0;
   keep transcript_id flag_cd4_75perc_apn5 flag_cd8_75perc_apn5 flag_cd19_75perc_apn5;
run;

proc freq data=flag_xs_dd_Event noprint;
    tables flag_cd4_75perc_apn5*flag_cd8_75perc_apn5*flag_cd19_75perc_apn5 / out=dd_count;
proc print data=dd_count;
run;

/*
  flag_cd4_    flag_cd8_    flag_cd19_
   75perc_      75perc_       75perc_
     apn5         apn5         apn5       COUNT

      0            0             0        59870
      0            0             1         2250
      0            1             0          832
      0            1             1          577
      1            0             0          402
      1            0             1          231
      1            1             0         2071
      1            1             1        20120
*/


data event.t1d_flag_xscript_dd_75perc_apn5;
   set flag_xs_dd_Event;
   if sum(flag_cd4_75perc_apn5,flag_cd8_75perc_apn5,flag_cd19_75perc_apn5) = 0 then flag_xscript_dd=.;
   else if sum(flag_cd4_75perc_apn5,flag_cd8_75perc_apn5,flag_cd19_75perc_apn5) < 3 then flag_xscript_dd=1;
   else flag_xscript_dd=0;
run;

