
proc sort data=pb_nic_junc;
   by junc_index;
proc sort data=pb_nnc_junc;
   by junc_index;
proc sort data=star_novel;
   by junc_index;
proc sort data=star_unannot;
   by junc_index;
proc sort data=cat_unannot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data star2event_novel;
  merge star_novel (in=in1) star_unannot (in=in2) cat_unannot (in=in3) pb_nic_junc (in=in4)
        pb_nnc_junc (in=in5) detection2;
  by junc_index;
  if in4 or in5 then flag_pb_junc=1; else flag_pb_junc=0;
  if in1 or in2 then flag_star_junc=1; else flag_star_junc=0;
  if in3 then flag_cat_junc=1; else flag_cat_junc=0;
  if in1 then flag_star_nnc=1; else flag_star_nnc=0;
  if in1 or in2 or in3 then output;
run;

* Check overlap;

*no NNC;
proc freq data=star2event_novel;
   where flag_star_nnc=0;
   tables flag_cat_junc*flag_star_junc;
run;

/*

 flag_cat_junc     flag_star_junc

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |   2311 |   2311
          |   0.00 |   4.20 |   4.20
          |   0.00 | 100.00 |
          |   0.00 |  29.30 |
 ---------+--------+--------+
        1 |  47142 |   5576 |  52718
          |  85.67 |  10.13 |  95.80
          |  89.42 |  10.58 |
          | 100.00 |  70.70 |
 ---------+--------+--------+
 Total       47142     7887    55029
             85.67    14.33   100.00


5576/7887 (70.70%) unannotated STAR junctions are also unannotated catalog junctions
5576/52718 (10.58%) unannotated catalog junctions are also in STAR

*/

*include NNC;
proc freq data=star2event_novel;
   tables flag_cat_junc*flag_star_junc;
run;


/*
  flag_cat_junc     flag_star_junc

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |  38875 |  38875
           |   0.00 |  42.44 |  42.44
           |   0.00 | 100.00 |
           |   0.00 |  87.46 |
  ---------+--------+--------+
         1 |  47142 |   5576 |  52718
           |  51.47 |   6.09 |  57.56
           |  89.42 |  10.58 |
           | 100.00 |  12.54 |
  ---------+--------+--------+
  Total       47142    44451    91593
              51.47    48.53   100.00


5576/44451 (12.54%) unannotated STAR junctions are also unannotated catalog junctions
5576/52718 (10.58%) unannotated catalog junctions are also in STAR


*/

/* NIC only */

* All;
proc freq data=star2event_novel;
   where flag_star_nnc=0;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*
flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  36427 |  16353 |  52780
         |  66.20 |  29.72 |  95.91
         |  69.02 |  30.98 |
         |  97.39 |  92.79 |
---------+--------+--------+
       1 |    978 |   1271 |   2249
         |   1.78 |   2.31 |   4.09
         |  43.49 |  56.51 |
         |   2.61 |   7.21 |
---------+--------+--------+
Total       37405    17624    55029
            67.97    32.03   100.00

TP = In both    = 1271
TN = In neither = 36427
FP = EA only    = 16353
FN = STAR only  = 978

FPR=FP/(FP+TN)=16353/(16353+36427)= 30.98%
FNR=FN/(FN+TP)=978/(978+1271)= 43.49%

flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  52900 |   1609 |  54509
         |  96.13 |   2.92 |  99.06
         |  97.05 |   2.95 |
         |  99.21 |  94.26 |
---------+--------+--------+
       1 |    422 |     98 |    520
         |   0.77 |   0.18 |   0.94
         |  81.15 |  18.85 |
         |   0.79 |   5.74 |
---------+--------+--------+
Total       53322     1707    55029
            96.90     3.10   100.00


TP = In both    = 98
TN = In neither = 52900
FP = EA only    = 1609
FN = STAR only  = 422

FPR=FP/(FP+TN)=1609/(1609+52900)=2.95%
FNR=FN/(FN+TP)=1079/(422+98)=81.15%


*/

* PB only;
proc freq data=star2event_novel;
   where flag_pb_junc=1 and flag_star_nnc=0;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*

flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |    136 |     24 |    160
         |  31.19 |   5.50 |  36.70
         |  85.00 |  15.00 |
         |  57.38 |  12.06 |
---------+--------+--------+
       1 |    101 |    175 |    276
         |  23.17 |  40.14 |  63.30
         |  36.59 |  63.41 |
         |  42.62 |  87.94 |
---------+--------+--------+
Total         237      199      436
            54.36    45.64   100.00

TP = In both    = 175
TN = In neither = 136
FP = EA only    = 24
FN = STAR only  = 101

FPR=FP/(FP+TN)=24/(24+136)= 15.00%
FNR=FN/(FN+TP)=101/(101+175)= 36.59%

flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |    306 |      1 |    307
         |  70.18 |   0.23 |  70.41
         |  99.67 |   0.33 |
         |  75.56 |   3.23 |
---------+--------+--------+
       1 |     99 |     30 |    129
         |  22.71 |   6.88 |  29.59
         |  76.74 |  23.26 |
         |  24.44 |  96.77 |
---------+--------+--------+
Total         405       31      436
            92.89     7.11   100.00



TP = In both    = 30
TN = In neither = 306
FP = EA only    = 1
FN = STAR only  = 99

FPR=FP/(FP+TN)=1/(1+306)= 0.33%
FNR=FN/(FN+TP)=99/(99+30)= 76.74%

*/

* Gene has PB;
proc freq data=star2event_novel;
   where flag_gene_has_pb=1 and flag_star_nnc=0;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;

/*

 flag_NPC_both_depth0
           flag_NPC_both_apn0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  23187 |  10789 |  33976
          |  64.96 |  30.22 |  95.18
          |  68.25 |  31.75 |
          |  96.88 |  91.72 |
 ---------+--------+--------+
        1 |    746 |    974 |   1720
          |   2.09 |   2.73 |   4.82
          |  43.37 |  56.63 |
          |   3.12 |   8.28 |
 ---------+--------+--------+
 Total       23933    11763    35696
             67.05    32.95   100.00


TP = In both    = 974
TN = In neither = 23187
FP = EA only    = 10789
FN = STAR only  = 746

FPR=FP/(FP+TN)=10789/(10789+23187)= 31.75%
FNR=FN/(FN+TP)=746/(746+974)= 43.37%

flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  34509 |    787 |  35296
         |  96.67 |   2.20 |  98.88
         |  97.77 |   2.23 |
         |  99.03 |  92.59 |
---------+--------+--------+
       1 |    337 |     63 |    400
         |   0.94 |   0.18 |   1.12
         |  84.25 |  15.75 |
         |   0.97 |   7.41 |
---------+--------+--------+
Total       34846      850    35696
            97.62     2.38   100.00


TP = In both    = 63
TN = In neither = 34509
FP = EA only    = 787
FN = STAR only  = 337

FPR=FP/(FP+TN)=787/(787+34509)= 2.23%
FNR=FN/(FN+TP)=337/(337+63)= 84.25%
*/



/* Unannot EA vs NIC/NNC STAR */

* All;
proc freq data=star2event_novel;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*
flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  67589 |  16353 |  83942
         |  73.79 |  17.85 |  91.65
         |  80.52 |  19.48 |
         |  91.37 |  92.79 |
---------+--------+--------+
       1 |   6380 |   1271 |   7651
         |   6.97 |   1.39 |   8.35
         |  83.39 |  16.61 |
         |   8.63 |   7.21 |
---------+--------+--------+
Total       73969    17624    91593
            80.76    19.24   100.00

TP = In both    = 1271
TN = In neither = 67589
FP = EA only    = 16353
FN = STAR only  = 6380

FPR=FP/(FP+TN)=16353/(16353+2734891)= 19.48%
FNR=FN/(FN+TP)=6380/(6380+1271)= 83.39%

flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  88807 |   1609 |  90416
         |  96.96 |   1.76 |  98.71
         |  98.22 |   1.78 |
         |  98.80 |  94.26 |
---------+--------+--------+
       1 |   1079 |     98 |   1177
         |   1.18 |   0.11 |   1.29
         |  91.67 |   8.33 |
         |   1.20 |   5.74 |
---------+--------+--------+
Total       89886     1707    91593
            98.14     1.86   100.00

TP = In both    = 98
TN = In neither = 88807
FP = EA only    = 1609
FN = STAR only  = 1079

FPR=FP/(FP+TN)=1609/(1609+88807)=1.78%
FNR=FN/(FN+TP)=1079/(1079+98)=91.67%


*/

* PB only;
proc freq data=star2event_novel;
   where flag_pb_junc=1;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*

flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |    530 |     24 |    554
         |  40.27 |   1.82 |  42.10
         |  95.67 |   4.33 |
         |  47.45 |  12.06 |
---------+--------+--------+
       1 |    587 |    175 |    762
         |  44.60 |  13.30 |  57.90
         |  77.03 |  22.97 |
         |  52.55 |  87.94 |
---------+--------+--------+
Total        1117      199     1316
            84.88    15.12   100.00

TP = In both    = 175
TN = In neither = 530
FP = EA only    = 24
FN = STAR only  = 587

FPR=FP/(FP+TN)=24/(24+530)= 4.33%
FNR=FN/(FN+TP)=587/(587+175)= 77.03%



flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   1039 |      1 |   1040
         |  78.95 |   0.08 |  79.03
         |  99.90 |   0.10 |
         |  80.86 |   3.23 |
---------+--------+--------+
       1 |    246 |     30 |    276
         |  18.69 |   2.28 |  20.97
         |  89.13 |  10.87 |
         |  19.14 |  96.77 |
---------+--------+--------+
Total        1285       31     1316
            97.64     2.36   100.00


TP = In both    = 30
TN = In neither = 1039
FP = EA only    = 1
FN = STAR only  = 246

FPR=FP/(FP+TN)=1/(1+1311)= 0.10%
FNR=FN/(FN+TP)=246/(246+30)= 89.13%

*/

* Gene has PB;
proc freq data=star2event_novel;
   where flag_gene_has_pb=1;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;

/*

flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  38108 |  10789 |  48897
         |  70.82 |  20.05 |  90.87
         |  77.94 |  22.06 |
         |  90.63 |  91.72 |
---------+--------+--------+
       1 |   3939 |    974 |   4913
         |   7.32 |   1.81 |   9.13
         |  80.18 |  19.82 |
         |   9.37 |   8.28 |
---------+--------+--------+
Total       42047    11763    53810
            78.14    21.86   100.00

TP = In both    = 974
TN = In neither = 38108
FP = EA only    = 10789
FN = STAR only  = 3939

FPR=FP/(FP+TN)=10789/(10789+38108)= 22.06%
FNR=FN/(FN+TP)=3939/(3939+974)= 80.18%


flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  52262 |    787 |  53049
         |  97.12 |   1.46 |  98.59
         |  98.52 |   1.48 |
         |  98.68 |  92.59 |
---------+--------+--------+
       1 |    698 |     63 |    761
         |   1.30 |   0.12 |   1.41
         |  91.72 |   8.28 |
         |   1.32 |   7.41 |
---------+--------+--------+
Total       52960      850    53810
            98.42     1.58   100.00


TP = In both    = 63
TN = In neither = 52262
FP = EA only    = 787
FN = STAR only  = 698

FPR=FP/(FP+TN)=787/(787+52262)= 1.48%
FNR=FN/(FN+TP)=698/(698+63)= 91.72%
*/


