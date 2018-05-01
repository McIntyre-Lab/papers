
proc sort data=pb_annot_junc;
   by junc_index;
proc sort data=star_annot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data pb2star_annot;
  merge pb_annot_junc (in=in1) star_annot (in=in2) detection2 (in=in3);
  by junc_index;
  if in1 then flag_junc_in_pacbio=1; else flag_junc_in_pacbio=0;
  if in2 then flag_junc_in_star=1; else flag_junc_in_star=0;
  if in1 or in2;
run;


proc freq data=pb2star_annot;
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
       0 |      0 |  69255 |  69255
         |   0.00 |  50.31 |  50.31
         |   0.00 | 100.00 |
         |   0.00 |  51.67 |
---------+--------+--------+
       1 |   3626 |  64767 |  68393
         |   2.63 |  47.05 |  49.69
         |   5.30 |  94.70 |
         | 100.00 |  48.33 |
---------+--------+--------+
Total        3626   134022   137648
             2.63    97.37   100.00


48.33% of NIC STAR junctions are also in PacBio
5.30% of NIC PacBio junctions are not in STAR set

*/

* all;
proc freq data=pb2star_annot;
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
       0 |  27324 |  41931 |  69255
         |  19.85 |  30.46 |  50.31
         |  39.45 |  60.55 |
         |  80.44 |  40.44 |
---------+--------+--------+
       1 |   6643 |  61750 |  68393
         |   4.83 |  44.86 |  49.69
         |   9.71 |  90.29 |
         |  19.56 |  59.56 |
---------+--------+--------+
Total       33967   103681   137648
            24.68    75.32   100.00




TP= In both          = 61750
TN= In neither       = 27324
FP= In STAR only     = 41931
FN= In PacBio only   = 6643

FPR=FP/(FP+TN)= 41931/(41931+27324)= 60.55%
FNR=FN/(FN+TP)= 6643/(6643+61750)= 9.71%



  flag_junc_in_pacbio
            flag_NPC_both_depth5

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  54379 |  14876 |  69255
           |  39.51 |  10.81 |  50.31
           |  78.52 |  21.48 |
           |  76.79 |  22.26 |
  ---------+--------+--------+
         1 |  16439 |  51954 |  68393
           |  11.94 |  37.74 |  49.69
           |  24.04 |  75.96 |
           |  23.21 |  77.74 |
  ---------+--------+--------+
  Total       70818    66830   137648
              51.45    48.55   100.00


TP= In both          = 51954
TN= In neither       = 54379
FP= In STAR only     = 14876
FN= In PacBio only   = 16439

FPR=FP/(FP+TN)= 14876/(14876+54379)= 21.48%
FNR=FN/(FN+TP)= 16439/(16439+51954)= 24.04%
*/

* PB junc only;
proc freq data=pb2star_annot;
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
proc freq data=pb2star_annot;
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
       0 |   5690 |   9999 |  15689
         |   6.77 |  11.89 |  18.66
         |  36.27 |  63.73 |
         |  46.14 |  13.94 |
---------+--------+--------+
       1 |   6643 |  61750 |  68393
         |   7.90 |  73.44 |  81.34
         |   9.71 |  90.29 |
         |  53.86 |  86.06 |
---------+--------+--------+
Total       12333    71749    84082
            14.67    85.33   100.00

           The SAS System       08:20 Thurs



TP= In both          = 61750
TN= In neither       = 5690
FP= In STAR only     = 9999
FN= In PacBio only   = 6643


FPR=FP/(FP+TN)=9999/(9999+5690)=63.73%
FNR=FN/(FN+TP)=6643/(6643+61750)=9.71%

    flag_junc_in_pacbio
              flag_NPC_both_depth5

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  12496 |   3193 |  15689
             |  14.86 |   3.80 |  18.66
             |  79.65 |  20.35 |
             |  43.19 |   5.79 |
    ---------+--------+--------+
           1 |  16439 |  51954 |  68393
             |  19.55 |  61.79 |  81.34
             |  24.04 |  75.96 |
             |  56.81 |  94.21 |
    ---------+--------+--------+
    Total       28935    55147    84082
                34.41    65.59   100.00


TP= In both          = 51954
TN= In neither       = 12496
FP= In STAR only     = 3193
FN= In PacBio only   = 16439


FPR=FP/(FP+TN)=3193/(3193+12496)=20.35%
FNR=FN/(FN+TP)=16439/(16439+51954)=24.04%
*/

