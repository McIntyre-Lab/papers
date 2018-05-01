
proc sort data=pb_annot_junc;
   by junc_index;
proc sort data=star_annot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data pb2star_nic;
  merge pb_annot_junc (in=in1) star_annot (in=in2) detection2 (in=in3);
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
       0 |      0 |  71808 |  71808
         |   0.00 |  51.22 |  51.22
         |   0.00 | 100.00 |
         |   0.00 |  52.06 |
---------+--------+--------+
       1 |   2276 |  66117 |  68393
         |   1.62 |  47.16 |  48.78
         |   3.33 |  96.67 |
         | 100.00 |  47.94 |
---------+--------+--------+
Total        2276   137925   140201
             1.62    98.38   100.00




47.94% of NIC STAR junctions are also in PacBio
3.33% of NIC PacBio junctions are not in STAR set

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
       0 |  29877 |  41931 |  71808
         |  21.31 |  29.91 |  51.22
         |  41.61 |  58.39 |
         |  81.81 |  40.44 |
---------+--------+--------+
       1 |   6643 |  61750 |  68393
         |   4.74 |  44.04 |  48.78
         |   9.71 |  90.29 |
         |  18.19 |  59.56 |
---------+--------+--------+
Total       36520   103681   140201
            26.05    73.95   100.00

           The SAS System       08:



TP= In both          = 61750
TN= In neither       = 29877
FP= In STAR only     = 41931
FN= In PacBio only   = 6643

FPR=FP/(FP+TN)= 41931/(41931+29877)= 58.39%
FNR=FN/(FN+TP)= 6643/(6643+61750)= 9.71%



 flag_junc_in_pacbio
           flag_NPC_both_depth5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  56932 |  14876 |  71808
          |  40.61 |  10.61 |  51.22
          |  79.28 |  20.72 |
          |  77.59 |  22.26 |
 ---------+--------+--------+
        1 |  16439 |  51954 |  68393
          |  11.73 |  37.06 |  48.78
          |  24.04 |  75.96 |
          |  22.41 |  77.74 |
 ---------+--------+--------+
 Total       73371    66830   140201
             52.33    47.67   100.00

TP= In both          = 51954
TN= In neither       = 56932
FP= In STAR only     = 14876
FN= In PacBio only   = 16439

FPR=FP/(FP+TN)= 14876/(14876+56932)= 20.72%
FNR=FN/(FN+TP)= 16439/(16439+51954)= 24.04%
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
       1 |   6643 |  61750 |  68393
         |   9.71 |  90.29 | 100.00
         |   9.71 |  90.29 |
         | 100.00 | 100.00 |
---------+--------+--------+
Total        6643    61750    68393
             9.71    90.29   100.00


TP= In both          = 61750
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 6643


FPR=FP/(FP+TN)= n/a
FNR=FN/(FN+TP)= 6643/(6643+61750)= 9.71%

   flag_junc_in_pacbio
             flag_NPC_both_depth5

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          1 |  16439 |  51954 |  68393
            |  24.04 |  75.96 | 100.00
            |  24.04 |  75.96 |
            | 100.00 | 100.00 |
   ---------+--------+--------+
   Total       16439    51954    68393
               24.04    75.96   100.00

TP= In both          = 51954
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 16439

FPR=FP/(FP+TN)=n/a
FNR=FN/(FN+TP)= 16439/(16439+51954)= 24.04%
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
       0 |   6083 |   9999 |  16082
         |   7.20 |  11.84 |  19.04
         |  37.82 |  62.18 |
         |  47.80 |  13.94 |
---------+--------+--------+
       1 |   6643 |  61750 |  68393
         |   7.86 |  73.10 |  80.96
         |   9.71 |  90.29 |
         |  52.20 |  86.06 |
---------+--------+--------+
Total       12726    71749    84475
            15.06    84.94   100.00

           The SAS System       08:

TP= In both          = 61750
TN= In neither       = 6083
FP= In STAR only     = 9999
FN= In PacBio only   = 6643


FPR=FP/(FP+TN)=9999/(9999+6083)=62.18%
FNR=FN/(FN+TP)=6643/(6643+61750)=9.71%

flag_junc_in_pacbio
          flag_NPC_both_depth5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  12889 |   3193 |  16082
         |  15.26 |   3.78 |  19.04
         |  80.15 |  19.85 |
         |  43.95 |   5.79 |
---------+--------+--------+
       1 |  16439 |  51954 |  68393
         |  19.46 |  61.50 |  80.96
         |  24.04 |  75.96 |
         |  56.05 |  94.21 |
---------+--------+--------+
Total       29328    55147    84475
            34.72    65.28   100.00




TP= In both          = 51954
TN= In neither       = 12889
FP= In STAR only     = 3193
FN= In PacBio only   = 16439


FPR=FP/(FP+TN)=3193/(3193+12889)=19.85%
FNR=FN/(FN+TP)=16439/(16439+51954)=24.04%
*/

