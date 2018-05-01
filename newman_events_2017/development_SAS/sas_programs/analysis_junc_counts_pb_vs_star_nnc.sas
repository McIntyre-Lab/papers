
proc sort data=pb_nnc_junc;
   by junc_index;
proc sort data=star_novel;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data pb2star_nnc;
  merge pb_nnc_junc (in=in1) star_novel (in=in2) detection2 (in=in3);
  by junc_index;
  if in1 then flag_junc_in_pacbio=1; else flag_junc_in_pacbio=0;
  if in2 then flag_junc_in_star=1; else flag_junc_in_star=0;
  if in1 or in2;
run;


proc freq data=pb2star_nnc;
  tables flaG_junc_in_pacbio*flag_junc_in_star;
run;

/*
flag_junc_in_pacbio
          flag_junc_in_star

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |      0 |  49560 |  49560
         |   0.00 |  93.07 |  93.07
         |   0.00 | 100.00 |
         |   0.00 |  98.15 |
---------+--------+--------+
       1 |   2755 |    936 |   3691
         |   5.17 |   1.76 |   6.93
         |  74.64 |  25.36 |
         | 100.00 |   1.85 |
---------+--------+--------+
Total        2755    50496    53251
             5.17    94.83   100.00


1.85% of NNC STAR junctions are also in PacBio
74.64% of NNC PacBio junctions are not in STAR set

*/

* all;
proc freq data=pb2star_nnc;
    tables flag_junc_in_pacbio*flag_NPC_both_depth0
           flag_junc_in_pacbio*flag_NPC_both_depth5;
run;

/*
 flag_junc_in_pacbio
           flag_NPC_both_depth0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  44644 |   4916 |  49560
          |  83.84 |   9.23 |  93.07
          |  90.08 |   9.92 |
          |  93.30 |  91.00 |
 ---------+--------+--------+
        1 |   3205 |    486 |   3691
          |   6.02 |   0.91 |   6.93
          |  86.83 |  13.17 |
          |   6.70 |   9.00 |
 ---------+--------+--------+
 Total       47849     5402    53251
             89.86    10.14   100.00

TP= In both          = 486
TN= In neither       = 44644
FP= In STAR only     = 4916
FN= In PacBio only   = 3205

FPR=FP/(FP+TN)= 4916/(4916+44644)= 9.92%
FNR=FN/(FN+TP)= 3205/(3205+486)= 86.83%

flag_junc_in_pacbio
          flag_NPC_both_depth5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  49050 |    510 |  49560
         |  92.11 |   0.96 |  93.07
         |  98.97 |   1.03 |
         |  93.26 |  77.63 |
---------+--------+--------+
       1 |   3544 |    147 |   3691
         |   6.66 |   0.28 |   6.93
         |  96.02 |   3.98 |
         |   6.74 |  22.37 |
---------+--------+--------+
Total       52594      657    53251
            98.77     1.23   100.00

TP= In both          = 147
TN= In neither       = 49050
FP= In STAR only     = 510
FN= In PacBio only   = 3544

FPR=FP/(FP+TN)= 510/(510+49050)= 1.03%
FNR=FN/(FN+TP)= 3544/(3544+147)= 96.02%
*/

* PB junc only;
proc freq data=pb2star_nnc;
    where flag_junc_in_pacbio=1;
    tables flag_junc_in_pacbio*flag_NPC_both_depth0
           flag_junc_in_pacbio*flag_NPC_both_depth5;
run;

/*
flag_junc_in_pacbio
          flag_NPC_both_depth0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       1 |   3205 |    486 |   3691
         |  86.83 |  13.17 | 100.00
         |  86.83 |  13.17 |
         | 100.00 | 100.00 |
---------+--------+--------+
Total        3205      486     3691
            86.83    13.17   100.00

TP= In both          = 486
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 3205

FPR=FP/(FP+TN)= n/a
FNR=FN/(FN+TP)= 3205/(3205+486)= 86.83%

 flag_junc_in_pacbio
           flag_NPC_both_depth5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        1 |   3544 |    147 |   3691
          |  96.02 |   3.98 | 100.00
          |  96.02 |   3.98 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total        3544      147     3691
             96.02     3.98   100.00

TP= In both          = 147
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 3544

FPR=FP/(FP+TN)=n/a
FNR=FN/(FN+TP)= 3544/(3544+147)= 96.02%
*/

* Gene has PB;
proc freq data=pb2star_nnc;
    where flag_gene_has_pb=1;
    tables flag_junc_in_pacbio*flag_NPC_both_depth0
           flag_junc_in_pacbio*flag_NPC_both_depth5;
run;

/*
 flag_junc_in_pacbio
           flag_NPC_both_depth0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  16468 |   2707 |  19175
          |  72.02 |  11.84 |  83.86
          |  85.88 |  14.12 |
          |  83.71 |  84.78 |
 ---------+--------+--------+
        1 |   3205 |    486 |   3691
          |  14.02 |   2.13 |  16.14
          |  86.83 |  13.17 |
          |  16.29 |  15.22 |
 ---------+--------+--------+
 Total       19673     3193    22866
             86.04    13.96   100.00


TP= In both          = 486
TN= In neither       = 16468
FP= In STAR only     = 2707
FN= In PacBio only   = 3205


FPR=FP/(FP+TN)=2707/(2707+16468)=14.12%
FNR=FN/(FN+TP)=3205/(3205+486)=86.83%

flag_junc_in_pacbio
          flag_NPC_both_depth5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  18961 |    214 |  19175
         |  82.92 |   0.94 |  83.86
         |  98.88 |   1.12 |
         |  84.25 |  59.28 |
---------+--------+--------+
       1 |   3544 |    147 |   3691
         |  15.50 |   0.64 |  16.14
         |  96.02 |   3.98 |
         |  15.75 |  40.72 |
---------+--------+--------+
Total       22505      361    22866
            98.42     1.58   100.00


TP= In both          = 147
TN= In neither       = 18961
FP= In STAR only     = 214
FN= In PacBio only   = 3544


FPR=FP/(FP+TN)=214/(214+18961)=1.12%
FNR=FN/(FN+TP)=3544/(3544+147)=96.02%
*/


