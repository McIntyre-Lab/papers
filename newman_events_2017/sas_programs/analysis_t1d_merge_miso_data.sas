
ods listing; ods html close;
libname con '!PATCON/sas_data';
libname event '!MCLAB/event_analysis/sas_data';

/* Merge event2gene with individual cell type MISO summaries and flag which events are seen in what samples */

data event2gene;
   set event.hg19_miso_event2gene;
   rename event_id=event_name;
run;

data cd4;
  set event.t1d_miso_summary_cd4;
  keep event_name;
run;

data cd8;
  set event.t1d_miso_summary_cd8;
  keep event_name;
run;

data cd19;
  set event.t1d_miso_summary_cd19;
  keep event_name;
run;

proc sort data=event2gene;
  by event_name;
proc sort data=cd4;
  by event_name;
proc sort data=cd8;
  by event_name;
proc sort data=cd19;
  by event_name;
run;

data flag_miso_event_dtct;
   merge event2gene (in=in1) cd4 (in=in2) cd8 (in=in3) cd19 (in=in4);
   by event_name;
   if in2 then flag_cd4_on=1; else flag_cd4_on=0;
   if in3 then flag_cd8_on=1; else flag_cd8_on=0;
   if in4 then flag_cd19_on=1; else flag_cd19_on=0;
   if in1 then output;
run;

proc freq data=flag_miso_event_dtct noprint;
  tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=check;
proc print data=check;
run;

/*

   flag_     flag_     flag_
  cd4_on    cd8_on    cd19_on    COUNT

     0         0         0       19176
     0         0         1        1326
     0         1         0         201
     0         1         1         151
     1         0         0        2902
     1         0         1        2694
     1         1         0        1231
     1         1         1       11541

*/


/* For each comparison, flag if diff >= 0.2 and BF>=5 or BF>=10, then merge with detection list */

data cd4cd8;
   set event.t1d_miso_compare_cd4cd8;
   if abs(diff) ge 0.2 then flag_cd4cd8_diff_ge_02=1; else flag_cd4cd8_diff_ge_02=0;
   if abs(bayes_factor) ge 10 then flag_cd4cd8_bf10=1; else flag_cd4cd8_bf10=0;
   if abs(bayes_factor) ge 5 then flag_cd4cd8_bf5=1; else flag_cd4cd8_bf5=0;
   keep event_name flag_cd4cd8_diff_ge_02 flag_cd4cd8_bf5 flag_cd4cd8_bf10;
run;

data cd4cd19;
   set event.t1d_miso_compare_cd4cd19;
   if abs(diff) ge 0.2 then flag_cd4cd19_diff_ge_02=1; else flag_cd4cd19_diff_ge_02=0;
   if abs(bayes_factor) ge 10 then flag_cd4cd19_bf10=1; else flag_cd4cd19_bf10=0;
   if abs(bayes_factor) ge 5 then flag_cd4cd19_bf5=1; else flag_cd4cd19_bf5=0;
   keep event_name flag_cd4cd19_diff_ge_02 flag_cd4cd19_bf5 flag_cd4cd19_bf10;
run;

data cd8cd19;
   set event.t1d_miso_compare_cd8cd19;
   if abs(diff) ge 0.2 then flag_cd8cd19_diff_ge_02=1; else flag_cd8cd19_diff_ge_02=0;
   if abs(bayes_factor) ge 10 then flag_cd8cd19_bf10=1; else flag_cd8cd19_bf10=0;
   if abs(bayes_factor) ge 5 then flag_cd8cd19_bf5=1; else flag_cd8cd19_bf5=0;
   keep event_name flag_cd8cd19_diff_ge_02 flag_cd8cd19_bf5 flag_cd8cd19_bf10;
run;

proc sort data=cd4cd8;
   by event_name;
proc sort data=cd4cd19;
   by event_name;
proc sort data=cd8cd19;
   by event_name;
proc sort data=flag_miso_event_dtct;
   by event_name;
run;

data t1d_miso_all;
   merge flag_miso_event_dtct (in=in1) cd4cd8 (in=in2)  cd4cd19 (in=in3)  cd8cd19 (in=in4) ;
   by event_name;
   if in2 then flag_cd4cd8_testable=1; else flag_cd4cd8_testable=0;
   if in3 then flag_cd4cd19_testable=1; else flag_cd4cd19_testable=0;
   if in4 then flag_cd8cd19_testable=1; else flag_cd8cd19_testable=0;
   if in1 then output;
run;

proc freq noprint data=t1d_miso_all;
   tables flag_cd4_on*flag_Cd8_on*flag_Cd19_on*
          flag_Cd4cd8_testable*flag_cd4cd19_testable*flag_cd8cd19_testable / out=test_check;
run;
proc print data=test_check;
run;

/*
       flag_    flag_    flag_    flag_cd4cd8_   flag_cd4cd19_   flag_cd8cd19_
Obs   cd4_on   cd8_on   cd19_on     testable        testable        testable     COUNT   PERCENT

 1       0        0        0            0              0               0         19176   48.8909
 2       0        0        1            0              0               0          1326    3.3808
 3       0        1        0            0              0               0           201    0.5125
 4       0        1        1            0              0               1           151    0.3850
 5       1        0        0            0              0               0          2902    7.3989

 6       1        0        1            0              1               0          2694    6.8686
 7       1        1        0            1              0               0          1231    3.1385
 8       1        1        1            1              1               1         11541   29.4248

*/


proc freq noprint data=t1d_miso_all;
    where flag_Cd4cd8_testable=1 or flag_cd4cd19_testable=1 or flag_cd8cd19_testable=1 ;

    tables flag_Cd4cd8_diff_ge_02*flag_Cd4cd8_bf5*flag_Cd4cd8_bf10*
          flag_cd4cd19_diff_ge_02*flag_cd4cd19_bf5*flag_cd4cd19_bf10*
          flag_cd8cd19_diff_ge_02*flag_cd8cd19_bf5*flag_cd8cd19_bf10 / out=test_check;
run;
proc print data=test_check;
run;

/* flag if different for any */

data t1d_miso_all2;
  set t1d_miso_all;
   if flag_cd4_on=1 or flag_cd8_on=1 or flag_Cd19_on=1 then flag_any_on=1; else flag_any_on=0;
   if sum(flag_cd4_on,flag_cd8_on,flag_Cd19_on)=1 then flag_1cell_event=1; else flag_1cell_event=0;
   if sum(flag_cd4_on,flag_cd8_on,flag_Cd19_on)=2 then flag_2cells_event=1; else flag_2cells_event=0;
   if flag_cd4_on=0 or flag_cd8_on=0 or flag_Cd19_on=0 then flag_all_off=1; else flag_all_off=0;
  if flag_cd4cd8_testable=1 or flag_cd4cd19_testable=1 or flag_cd8cd19_testable=1 then flag_miso_testable=1;
  else flag_miso_testable=0;
  if flag_Cd4cd8_diff_ge_02=1 and flag_Cd4cd8_bf5=1 then flag_sig_cd4cd8_bf5=1; else flag_sig_cd4cd8_bf5=0;
  if flag_Cd4cd8_diff_ge_02=1 and flag_cd4cd19_bf5=1 then flag_sig_cd4cd19_bf5=1; else flag_sig_cd4cd19_bf5=0;
  if flag_Cd4cd8_diff_ge_02=1 and flag_cd8cd19_bf5=1 then flag_sig_cd8cd19_bf5=1; else flag_sig_cd8cd19_bf5=0;
  if flag_sig_cd4cd8_bf5=1 or flag_sig_cd4cd8_bf5=1 or flag_sig_cd4cd8_bf5=1 then flag_sig_bf5=1;
  else flag_sig_bf5=0;

  if flag_Cd4cd8_diff_ge_02=1 and flag_Cd4cd8_bf10=1 then flag_sig_cd4cd8_bf10=1; else flag_sig_cd4cd8_bf10=0;
  if flag_Cd4cd8_diff_ge_02=1 and flag_cd4cd19_bf10=1 then flag_sig_cd4cd19_bf10=1; else flag_sig_cd4cd19_bf10=0;
  if flag_Cd4cd8_diff_ge_02=1 and flag_cd8cd19_bf10=1 then flag_sig_cd8cd19_bf10=1; else flag_sig_cd8cd19_bf10=0;
  if flag_sig_cd4cd8_bf10=1 or flag_sig_cd4cd8_bf10=1 or flag_sig_cd4cd8_bf10=1 then flag_sig_bf10=1;
  else flag_sig_bf10=0;
run;

proc freq data=t1d_miso_all2;
   tables flag_any_on*flag_miso_testable flag_sig_bf5*flag_sig_bf10;
run;


/*

 flag_any_on     flag_miso_testable

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  19176 |      0 |  19176
          |  48.89 |   0.00 |  48.89
          | 100.00 |   0.00 |
          |  81.24 |   0.00 |
 ---------+--------+--------+
        1 |   4429 |  15617 |  20046
          |  11.29 |  39.82 |  51.11
          |  22.09 |  77.91 |
          |  18.76 | 100.00 |
 ---------+--------+--------+
 Total       23605    15617    39222
             60.18    39.82   100.00


 flag_sig_bf5     flag_sig_bf10

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  39182 |      0 |  39182
          |  99.90 |   0.00 |  99.90
          | 100.00 |   0.00 |
          |  99.93 |   0.00 |
 ---------+--------+--------+
        1 |     27 |     13 |     40
          |   0.07 |   0.03 |   0.10
          |  67.50 |  32.50 |
          |   0.07 | 100.00 |
 ---------+--------+--------+
 Total       39209       13    39222
             99.97     0.03   100.00

*/

/* Make permenant */

data event.t1d_all_miso_results_v2;
  set t1d_miso_all2;
run;

