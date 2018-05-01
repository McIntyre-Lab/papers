
proc sort data=pb_Annot_junc;
   by junc_index;
proc sort data=star_Annot;
   by junc_index;
proc sort data=cat_Annot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data star2event_Annot;
  merge star_Annot (in=in1) cat_Annot (in=in2) pb_Annot_junc (in=in3) detection2;
  by junc_index;
  if in3 then flag_pb_junc=1; else flag_pb_junc=0;
  if in1 then flag_star_junc=1; else flag_star_junc=0;
  if in2 then flag_cat_junc=1; else flag_cat_junc=0;
  if in1 or in2 then output;
run;

* Check overlap;

proc freq data=star2event_Annot;
   tables flag_cat_junc*flag_star_junc;
run;

/*

Table of flag_cat_junc by flag_star_junc

  flag_cat_junc     flag_star_junc

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |  12742 |  12742
           |   0.00 |   9.12 |   9.12
           |   0.00 | 100.00 |
           |   0.00 |   9.51 |
  ---------+--------+--------+
         1 |   5655 | 121280 | 126935
           |   4.05 |  86.83 |  90.88
           |   4.46 |  95.54 |
           | 100.00 |  90.49 |
  ---------+--------+--------+
  Total        5655   134022   139677
               4.05    95.95   100.00


121280/134022 (90.49%) annotated STAR junctions are also annotated catalog junctions 
121280/126935 (95.54%) annotated catalog junctions are also in STAR

*/


proc freq data=star2event_Annot;
   where flag_pb_junc=1;
   tables flag_cat_junc*flag_star_junc;
run;


/*

  flag_cat_junc     flag_star_junc

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   2414 |   2414
           |   0.00 |   3.64 |   3.64
           |   0.00 | 100.00 |
           |   0.00 |   3.73 |
  ---------+--------+--------+
         1 |   1567 |  62353 |  63920
           |   2.36 |  94.00 |  96.36
           |   2.45 |  97.55 |
           | 100.00 |  96.27 |
  ---------+--------+--------+
  Total        1567    64767    66334
               2.36    97.64   100.00
*/

* All;
proc freq data=star2event_annot;
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
       0 |  31271 |   4725 |  35996
         |  22.39 |   3.38 |  25.77
         |  86.87 |  13.13 |
         |  70.60 |   4.95 |
---------+--------+--------+
       1 |  13024 |  90657 | 103681
         |   9.32 |  64.90 |  74.23
         |  12.56 |  87.44 |
         |  29.40 |  95.05 |
---------+--------+--------+
Total       44295    95382   139677
            31.71    68.29   100.00



TP = In both    = 90657
TN = In neither = 31271
FP = EA only    = 4725
FN = STAR only  = 13024

FPR=FP/(FP+TN)=4725/(4725+31271)= 13.13%
FNR=FN/(FN+TP)=13024/(13024+90657)= 12.56%

 flag_NPC_both_depth5
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  70761 |   2086 |  72847
          |  50.66 |   1.49 |  52.15
          |  97.14 |   2.86 |
          |  75.72 |   4.51 |
 ---------+--------+--------+
        1 |  22685 |  44145 |  66830
          |  16.24 |  31.61 |  47.85
          |  33.94 |  66.06 |
          |  24.28 |  95.49 |
 ---------+--------+--------+
 Total       93446    46231   139677
             66.90    33.10   100.00


TP = In both    = 44145
TN = In neither = 70761
FP = EA only    = 2086
FN = STAR only  = 22685

FPR=FP/(FP+TN)=2086/(2086+70761)=2.86%
FNR=FN/(FN+TP)=22685/(22685+44145)=33.94%


*/

* PB only;
proc freq data=star2event_annot;
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
        0 |   2600 |   1984 |   4584
          |   3.92 |   2.99 |   6.91
          |  56.72 |  43.28 |
          |  41.73 |   3.30 |
 ---------+--------+--------+
        1 |   3631 |  58119 |  61750
          |   5.47 |  87.62 |  93.09
          |   5.88 |  94.12 |
          |  58.27 |  96.70 |
 ---------+--------+--------+
 Total        6231    60103    66334
              9.39    90.61   100.00

TP = In both    = 58119
TN = In neither = 2600
FP = EA only    = 1984
FN = STAR only  = 3631

FPR=FP/(FP+TN)=1984/(1984+2600)= 43.28%
FNR=FN/(FN+TP)=3631/(3631+58119)= 5.88%


 flag_NPC_both_depth5
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  12904 |   1476 |  14380
          |  19.45 |   2.23 |  21.68
          |  89.74 |  10.26 |
          |  48.45 |   3.72 |
 ---------+--------+--------+
        1 |  13727 |  38227 |  51954
          |  20.69 |  57.63 |  78.32
          |  26.42 |  73.58 |
          |  51.55 |  96.28 |
 ---------+--------+--------+
 Total       26631    39703    66334
             40.15    59.85   100.00


TP = In both    = 38227
TN = In neither = 12904
FP = EA only    = 1476
FN = STAR only  = 13727

FPR=FP/(FP+TN)=1476/(1476+12904)= 10.26%
FNR=FN/(FN+TP)=13727/(13727+38227)= 26.42%

*/

* Gene has PB;
proc freq data=star2event_annot;
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
       0 |   8673 |   2740 |  11413
         |  10.43 |   3.29 |  13.72
         |  75.99 |  24.01 |
         |  60.00 |   3.99 |
---------+--------+--------+
       1 |   5783 |  65966 |  71749
         |   6.95 |  79.32 |  86.28
         |   8.06 |  91.94 |
         |  40.00 |  96.01 |
---------+--------+--------+
Total       14456    68706    83162
            17.38    82.62   100.00

           The SAS System       08:20

TP = In both    = 65966
TN = In neither = 8673
FP = EA only    = 2740
FN = STAR only  = 5783

FPR=FP/(FP+TN)=2740/(2740+8673)= 24.01%
FNR=FN/(FN+TP)=5783/(5783+65966)= 8.06%


flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  26447 |   1568 |  28015
         |  31.80 |   1.89 |  33.69
         |  94.40 |   5.60 |
         |  62.54 |   3.84 |
---------+--------+--------+
       1 |  15841 |  39306 |  55147
         |  19.05 |  47.26 |  66.31
         |  28.73 |  71.27 |
         |  37.46 |  96.16 |
---------+--------+--------+
Total       42288    40874    83162
            50.85    49.15   100.00


TP = In both    = 39306
TN = In neither = 26447
FP = EA only    = 1568
FN = STAR only  = 15841

FPR=FP/(FP+TN)=1568/(1568+26447)= 5.60%
FNR=FN/(FN+TP)=15841/(15841+39306)= 28.73%
*/

