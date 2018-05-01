
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

 flag_cat_junc     flag_star_junc

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        1 | 152130 | 137925 | 290055
          |  52.45 |  47.55 | 100.00
          |  52.45 |  47.55 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total      152130   137925   290055
             52.45    47.55   100.00



137925/137925 (100%) annotated STAR junctions are also annotated catalog junctions (expected)
137925/290055 (47.55%) annotated catalog junctions are also in STAR

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
        0 | 181649 |   4725 | 186374
          |  62.63 |   1.63 |  64.25
          |  97.46 |   2.54 |
          |  93.31 |   4.95 |
 ---------+--------+--------+
        1 |  13024 |  90657 | 103681
          |   4.49 |  31.26 |  35.75
          |  12.56 |  87.44 |
          |   6.69 |  95.05 |
 ---------+--------+--------+
 Total      194673    95382   290055
             67.12    32.88   100.00

TP = In both    = 90657
TN = In neither = 181649
FP = EA only    = 4725
FN = STAR only  = 13024

FPR=FP/(FP+TN)=4725/(4725+181649)= 2.54%
FNR=FN/(FN+TP)=13024/(13024+90657)= 12.56%


 flag_NPC_both_depth5
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 221139 |   2086 | 223225
          |  76.24 |   0.72 |  76.96
          |  99.07 |   0.93 |
          |  90.70 |   4.51 |
 ---------+--------+--------+
        1 |  22685 |  44145 |  66830
          |   7.82 |  15.22 |  23.04
          |  33.94 |  66.06 |
          |   9.30 |  95.49 |
 ---------+--------+--------+
 Total      243824    46231   290055
             84.06    15.94   100.00




TP = In both    = 44145
TN = In neither = 221139
FP = EA only    = 2086
FN = STAR only  = 22685

FPR=FP/(FP+TN)=2086/(2086+221139)=0.93%
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
        0 |   4659 |   1984 |   6643
          |   6.81 |   2.90 |   9.71
          |  70.13 |  29.87 |
          |  56.20 |   3.30 |
 ---------+--------+--------+
        1 |   3631 |  58119 |  61750
          |   5.31 |  84.98 |  90.29
          |   5.88 |  94.12 |
          |  43.80 |  96.70 |
 ---------+--------+--------+
 Total        8290    60103    68393
             12.12    87.88   100.00


TP = In both    = 58119
TN = In neither = 4659
FP = EA only    = 1984
FN = STAR only  = 3631

FPR=FP/(FP+TN)=1984/(1984+4659)= 29.87%
FNR=FN/(FN+TP)=3631/(3631+58119)= 5.88%

  flag_NPC_both_depth5
            flag_NPC_both_apn5

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  14963 |   1476 |  16439
           |  21.88 |   2.16 |  24.04
           |  91.02 |   8.98 |
           |  52.15 |   3.72 |
  ---------+--------+--------+
         1 |  13727 |  38227 |  51954
           |  20.07 |  55.89 |  75.96
           |  26.42 |  73.58 |
           |  47.85 |  96.28 |
  ---------+--------+--------+
  Total       28690    39703    68393
              41.95    58.05   100.00


TP = In both    = 38227
TN = In neither = 14963
FP = EA only    = 1476
FN = STAR only  = 13727

FPR=FP/(FP+TN)=1476/(1476+14963)= 8.98%
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
       0 |  23161 |   2740 |  25901
         |  23.72 |   2.81 |  26.52
         |  89.42 |  10.58 |
         |  80.02 |   3.99 |
---------+--------+--------+
       1 |   5783 |  65966 |  71749
         |   5.92 |  67.55 |  73.48
         |   8.06 |  91.94 |
         |  19.98 |  96.01 |
---------+--------+--------+
Total       28944    68706    97650
            29.64    70.36   100.00

           The SAS System       08:20

TP = In both    = 65966
TN = In neither = 23161
FP = EA only    = 2740
FN = STAR only  = 5783

FPR=FP/(FP+TN)=2740/(2740+23161)= 10.58%
FNR=FN/(FN+TP)=5783/(5783+65966)= 8.06%

 flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  40935 |   1568 |  42503
         |  41.92 |   1.61 |  43.53
         |  96.31 |   3.69 |
         |  72.10 |   3.84 |
---------+--------+--------+
       1 |  15841 |  39306 |  55147
         |  16.22 |  40.25 |  56.47
         |  28.73 |  71.27 |
         |  27.90 |  96.16 |
---------+--------+--------+
Total       56776    40874    97650
            58.14    41.86   100.00


TP = In both    = 39306
TN = In neither = 40935
FP = EA only    = 1568
FN = STAR only  = 15841

FPR=FP/(FP+TN)=1568/(1568+40935)= 3.69%
FNR=FN/(FN+TP)=15841/(15841+39306)= 28.73%
*/



