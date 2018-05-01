
/* Compare STAR to PacBio -- unannotated junctions */
proc sort data=pb_nic_junc;
   by junc_index;
proc sort data=star_unannot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data pb2star_nic;
  merge pb_nic_junc (in=in1) star_unannot (in=in2) detection2 (in=in3);
  by junc_index;
  if in1 then flag_junc_in_pacbio=1; else flag_junc_in_pacbio=0;
  if in2 then flag_junc_in_star=1; else flag_junc_in_star=0;
  if in1 or in2;
run;



proc freq data=pb2star_nic;
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
       0 |      0 |   9151 |   9151
         |   0.00 |  93.35 |  93.35
         |   0.00 | 100.00 |
         |   0.00 |  95.68 |
---------+--------+--------+
       1 |    239 |    413 |    652
         |   2.44 |   4.21 |   6.65
         |  36.66 |  63.34 |
         | 100.00 |   4.32 |
---------+--------+--------+
Total         239     9564     9803
             2.44    97.56   100.00

4.32% of NIC STAR junctions are also in PacBio
36.66% of NIC PacBio junctions are not in STAR set

*/

* all;
proc freq data=pb2star_nic;
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
        0 |   7178 |   1973 |   9151
          |  73.22 |  20.13 |  93.35
          |  78.44 |  21.56 |
          |  95.02 |  87.73 |
 ---------+--------+--------+
        1 |    376 |    276 |    652
          |   3.84 |   2.82 |   6.65
          |  57.67 |  42.33 |
          |   4.98 |  12.27 |
 ---------+--------+--------+
 Total        7554     2249     9803
             77.06    22.94   100.00



TP= In both          = 276
TN= In neither       = 7178
FP= In STAR only     = 1973
FN= In PacBio only   = 376

FPR=FP/(FP+TN)= 1973/(1973+7178)= 21.56%
FNR=FN/(FN+TP)= 376/(376+276)= 57.67%

 flag_junc_in_pacbio
           flag_NPC_both_depth5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   8760 |    391 |   9151
          |  89.36 |   3.99 |  93.35
          |  95.73 |   4.27 |
          |  94.37 |  75.19 |
 ---------+--------+--------+
        1 |    523 |    129 |    652
          |   5.34 |   1.32 |   6.65
          |  80.21 |  19.79 |
          |   5.63 |  24.81 |
 ---------+--------+--------+
 Total        9283      520     9803
             94.70     5.30   100.00

TP= In both          = 129
TN= In neither       = 8760
FP= In STAR only     = 391
FN= In PacBio only   = 523

FPR=FP/(FP+TN)= 391/(391+8760)= 4.27%
FNR=FN/(FN+TP)= 523/(523+129)= 80.21%
*/

* PB junc only;
proc freq data=pb2star_nic;
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
       1 |    376 |    276 |    652
         |  57.67 |  42.33 | 100.00
         |  57.67 |  42.33 |
         | 100.00 | 100.00 |
---------+--------+--------+
Total         376      276      652
            57.67    42.33   100.00


TP= In both          = 276
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 376


FPR=FP/(FP+TN)= n/a
FNR=FN/(FN+TP)= 376/(376+276)= 42.33%

 flag_junc_in_pacbio
           flag_NPC_both_depth5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        1 |    523 |    129 |    652
          |  80.21 |  19.79 | 100.00
          |  80.21 |  19.79 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total         523      129      652
             80.21    19.79   100.00

TP= In both          = 129
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 523

FPR=FP/(FP+TN)=n/a
FNR=FN/(FN+TP)= 523/(523+129)= 80.21%
*/

* Gene has PB;
proc freq data=pb2star_nic;
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
       0 |   5193 |   1444 |   6637
         |  71.24 |  19.81 |  91.06
         |  78.24 |  21.76 |
         |  93.25 |  83.95 |
---------+--------+--------+
       1 |    376 |    276 |    652
         |   5.16 |   3.79 |   8.94
         |  57.67 |  42.33 |
         |   6.75 |  16.05 |
---------+--------+--------+
Total        5569     1720     7289
            76.40    23.60   100.00

TP= In both          = 276
TN= In neither       = 5193
FP= In STAR only     = 1444
FN= In PacBio only   = 376


FPR=FP/(FP+TN)=1444/(1444+5193)=21.76%
FNR=FN/(FN+TP)=376/(376+276)=57.67%


flag_junc_in_pacbio
          flag_NPC_both_depth5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   6366 |    271 |   6637
         |  87.34 |   3.72 |  91.06
         |  95.92 |   4.08 |
         |  92.41 |  67.75 |
---------+--------+--------+
       1 |    523 |    129 |    652
         |   7.18 |   1.77 |   8.94
         |  80.21 |  19.79 |
         |   7.59 |  32.25 |
---------+--------+--------+
Total        6889      400     7289
            94.51     5.49   100.00

TP= In both          = 129
TN= In neither       = 6366
FP= In STAR only     = 271
FN= In PacBio only   = 523


FPR=FP/(FP+TN)=271/(271+6366)=4.08%
FNR=FN/(FN+TP)=523/(523+129)=80.21%
*/

