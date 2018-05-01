

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
        0 |      0 |   7483 |   7483
          |   0.00 |  91.99 |  91.99
          |   0.00 | 100.00 |
          |   0.00 |  94.88 |
 ---------+--------+--------+
        1 |    248 |    404 |    652
          |   3.05 |   4.97 |   8.01
          |  38.04 |  61.96 |
          | 100.00 |   5.12 |
 ---------+--------+--------+
 Total         248     7887     8135
              3.05    96.95   100.00

5.12% of NIC STAR junctions are also in PacBio
38.04% of NIC PacBio junctions are not in STAR set

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
        0 |   5510 |   1973 |   7483
          |  67.73 |  24.25 |  91.99
          |  73.63 |  26.37 |
          |  93.61 |  87.73 |
 ---------+--------+--------+
        1 |    376 |    276 |    652
          |   4.62 |   3.39 |   8.01
          |  57.67 |  42.33 |
          |   6.39 |  12.27 |
 ---------+--------+--------+
 Total        5886     2249     8135
             72.35    27.65   100.00

TP= In both          = 276
TN= In neither       = 5510
FP= In STAR only     = 1973
FN= In PacBio only   = 376

FPR=FP/(FP+TN)= 1973/(1973+5510)= 26.37%
FNR=FN/(FN+TP)= 376/(376+276)= 57.67%

 flag_junc_in_pacbio
           flag_NPC_both_depth5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   7092 |    391 |   7483
          |  87.18 |   4.81 |  91.99
          |  94.77 |   5.23 |
          |  93.13 |  75.19 |
 ---------+--------+--------+
        1 |    523 |    129 |    652
          |   6.43 |   1.59 |   8.01
          |  80.21 |  19.79 |
          |   6.87 |  24.81 |
 ---------+--------+--------+
 Total        7615      520     8135
             93.61     6.39   100.00


TP= In both          = 129
TN= In neither       = 7092
FP= In STAR only     = 391
FN= In PacBio only   = 523

FPR=FP/(FP+TN)= 391/(391+7092)= 5.23%
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
        0 |   4155 |   1444 |   5599
          |  66.47 |  23.10 |  89.57
          |  74.21 |  25.79 |
          |  91.70 |  83.95 |
 ---------+--------+--------+
        1 |    376 |    276 |    652
          |   6.02 |   4.42 |  10.43
          |  57.67 |  42.33 |
          |   8.30 |  16.05 |
 ---------+--------+--------+
 Total        4531     1720     6251
             72.48    27.52   100.00

TP= In both          = 276
TN= In neither       = 4155
FP= In STAR only     = 1444
FN= In PacBio only   = 376


FPR=FP/(FP+TN)=1444/(1444+5193)=25.79%
FNR=FN/(FN+TP)=376/(376+276)=57.67%

 flag_junc_in_pacbio
           flag_NPC_both_depth5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   5328 |    271 |   5599
          |  85.23 |   4.34 |  89.57
          |  95.16 |   4.84 |
          |  91.06 |  67.75 |
 ---------+--------+--------+
        1 |    523 |    129 |    652
          |   8.37 |   2.06 |  10.43
          |  80.21 |  19.79 |
          |   8.94 |  32.25 |
 ---------+--------+--------+
 Total        5851      400     6251
             93.60     6.40   100.00


TP= In both          = 129
TN= In neither       = 5328
FP= In STAR only     = 271
FN= In PacBio only   = 523


FPR=FP/(FP+TN)=271/(271+6366)=4.84%
FNR=FN/(FN+TP)=523/(523+129)=80.21%
*/

