
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
       0 |      0 |  35684 |  35684
         |   0.00 |  90.63 |  90.63
         |   0.00 | 100.00 |
         |   0.00 |  97.59 |
---------+--------+--------+
       1 |   2811 |    880 |   3691
         |   7.14 |   2.23 |   9.37
         |  76.16 |  23.84 |
         | 100.00 |   2.41 |
---------+--------+--------+
Total        2811    36564    39375
             7.14    92.86   100.00


2.41% of NNC STAR junctions are also in PacBio
76.16% of NNC PacBio junctions are not in STAR set

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
       0 |  30768 |   4916 |  35684
         |  78.14 |  12.49 |  90.63
         |  86.22 |  13.78 |
         |  90.57 |  91.00 |
---------+--------+--------+
       1 |   3205 |    486 |   3691
         |   8.14 |   1.23 |   9.37
         |  86.83 |  13.17 |
         |   9.43 |   9.00 |
---------+--------+--------+
Total       33973     5402    39375
            86.28    13.72   100.00

TP= In both          = 486
TN= In neither       = 30768
FP= In STAR only     = 4916
FN= In PacBio only   = 3205

FPR=FP/(FP+TN)= 4916/(4916+30768)= 13.78%
FNR=FN/(FN+TP)= 3205/(3205+486)= 86.83%



flag_junc_in_pacbio
          flag_NPC_both_depth5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  35174 |    510 |  35684
         |  89.33 |   1.30 |  90.63
         |  98.57 |   1.43 |
         |  90.85 |  77.63 |
---------+--------+--------+
       1 |   3544 |    147 |   3691
         |   9.00 |   0.37 |   9.37
         |  96.02 |   3.98 |
         |   9.15 |  22.37 |
---------+--------+--------+
Total       38718      657    39375
            98.33     1.67   100.00

TP= In both          = 147
TN= In neither       = 35174
FP= In STAR only     = 510
FN= In PacBio only   = 3544

FPR=FP/(FP+TN)= 510/(510+35174)= 1.43%
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
       0 |  14527 |   2707 |  17234
         |  69.42 |  12.94 |  82.36
         |  84.29 |  15.71 |
         |  81.93 |  84.78 |
---------+--------+--------+
       1 |   3205 |    486 |   3691
         |  15.32 |   2.32 |  17.64
         |  86.83 |  13.17 |
         |  18.07 |  15.22 |
---------+--------+--------+
Total       17732     3193    20925
            84.74    15.26   100.00

TP= In both          = 486
TN= In neither       = 14527
FP= In STAR only     = 2707
FN= In PacBio only   = 3205


FPR=FP/(FP+TN)=2707/(2707+14527)=15.71%
FNR=FN/(FN+TP)=3205/(3205+486)=86.83%

flag_junc_in_pacbio
          flag_NPC_both_depth5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  17020 |    214 |  17234
         |  81.34 |   1.02 |  82.36
         |  98.76 |   1.24 |
         |  82.77 |  59.28 |
---------+--------+--------+
       1 |   3544 |    147 |   3691
         |  16.94 |   0.70 |  17.64
         |  96.02 |   3.98 |
         |  17.23 |  40.72 |
---------+--------+--------+
Total       20564      361    20925
            98.27     1.73   100.00


TP= In both          = 147
TN= In neither       = 17020
FP= In STAR only     = 214
FN= In PacBio only   = 3544


FPR=FP/(FP+TN)=214/(214+17020)=1.24%
FNR=FN/(FN+TP)=3544/(3544+147)=96.02%
*/


