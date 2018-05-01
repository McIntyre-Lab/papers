
/* Compare STAR to PacBio -- annotated junctions */

proc sort data=pb_annot_junc;
   by junc_index;
proc sort data=cat_annot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data pb2cat_annot;
  merge pb_annot_junc (in=in1) cat_annot (in=in2) detection2 (in=in3);
  by junc_index;
  if in1 then flag_junc_in_pacbio=1; else flag_junc_in_pacbio=0;
  if in2 then flag_junc_in_catalog=1; else flag_junc_in_catalog=0;
  if in1 or in2;
run;


proc freq data=pb2cat_annot;
  tables flaG_junc_in_pacbio*flag_junc_in_catalog;
run;

/*
 flag_junc_in_pacbio
           flag_junc_in_catalog

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        0 | 221662 | 221662
          |  76.42 |  76.42
          | 100.00 |
          |  76.42 |
 ---------+--------+
        1 |  68393 |  68393
          |  23.58 |  23.58
          | 100.00 |
          |  23.58 |
 ---------+--------+
 Total      290055   290055
            100.00   100.00

100% of annotated PacBio junctions are in the junction catalog (expected)
Though, these only make up 23.58% of all annotated junctions

*/

* all;
proc freq data=pb2cat_annot;
    tables flag_junc_in_pacbio*flag_NPC_both_apn0
           flag_junc_in_pacbio*flag_NPC_both_apn5;
run;

/*

 flag_junc_in_pacbio
           flag_NPC_both_apn0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 186383 |  35279 | 221662
          |  64.26 |  12.16 |  76.42
          |  84.08 |  15.92 |
          |  95.74 |  36.99 |
 ---------+--------+--------+
        1 |   8290 |  60103 |  68393
          |   2.86 |  20.72 |  23.58
          |  12.12 |  87.88 |
          |   4.26 |  63.01 |
 ---------+--------+--------+
 Total      194673    95382   290055
             67.12    32.88   100.00

            The SAS System       08:20 T


TP= In both          = 60103
TN= In neither       = 186383
FP= In STAR only     = 35279
FN= In PacBio only   = 8290

FPR=FP/(FP+TN)= 35279/(35279+186383)= 15.92%
FNR=FN/(FN+TP)= 8290/(8290+60103)= 12.12%

 flag_junc_in_pacbio
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 215134 |   6528 | 221662
          |  74.17 |   2.25 |  76.42
          |  97.05 |   2.95 |
          |  88.23 |  14.12 |
 ---------+--------+--------+
        1 |  28690 |  39703 |  68393
          |   9.89 |  13.69 |  23.58
          |  41.95 |  58.05 |
          |  11.77 |  85.88 |
 ---------+--------+--------+
 Total      243824    46231   290055
             84.06    15.94   100.00



TP= In both          = 39703
TN= In neither       = 215134
FP= In STAR only     = 6528
FN= In PacBio only   = 28690

FPR=FP/(FP+TN)= 6528/(6528+215134)= 2.95%
FNR=FN/(FN+TP)= 28690/(28690+39703)= 41.95%
*/

* PB junc only;
proc freq data=pb2cat_annot;
    where flag_junc_in_pacbio=1;
    tables flag_junc_in_pacbio*flag_NPC_both_apn0
           flag_junc_in_pacbio*flag_NPC_both_apn5;
run;

/*

 flag_junc_in_pacbio
           flag_NPC_both_apn0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        1 |   8290 |  60103 |  68393
          |  12.12 |  87.88 | 100.00
          |  12.12 |  87.88 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total        8290    60103    68393
             12.12    87.88   100.00

TP= In both          = 60103
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 8290


FPR=FP/(FP+TN)= n/a
FNR=FN/(FN+TP)= 8290/(8290+60103)= 12.12%


flag_junc_in_pacbio
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       1 |  28690 |  39703 |  68393
         |  41.95 |  58.05 | 100.00
         |  41.95 |  58.05 |
         | 100.00 | 100.00 |
---------+--------+--------+
Total       28690    39703    68393
            41.95    58.05   100.00

TP= In both          = 39703
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 28690

FPR=FP/(FP+TN)=n/a
FNR=FN/(FN+TP)= 28690/(28690+39703)= 41.95%
*/

* Gene has PB;
proc freq data=pb2cat_annot;
    where flag_gene_has_pb=1;
    tables flag_junc_in_pacbio*flag_NPC_both_apn0
           flag_junc_in_pacbio*flag_NPC_both_apn5;
run;

/*
 flag_junc_in_pacbio
           flag_NPC_both_apn0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  20654 |   8603 |  29257
          |  21.15 |   8.81 |  29.96
          |  70.60 |  29.40 |
          |  71.36 |  12.52 |
 ---------+--------+--------+
        1 |   8290 |  60103 |  68393
          |   8.49 |  61.55 |  70.04
          |  12.12 |  87.88 |
          |  28.64 |  87.48 |
 ---------+--------+--------+
 Total       28944    68706    97650
             29.64    70.36   100.00

TP= In both          = 60103
TN= In neither       = 20654
FP= In STAR only     = 8603
FN= In PacBio only   = 8290


FPR=FP/(FP+TN)=8603/(8603+20654)=29.40%
FNR=FN/(FN+TP)=8290/(8290+60103)=12.12%


 flag_junc_in_pacbio
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  28086 |   1171 |  29257
          |  28.76 |   1.20 |  29.96
          |  96.00 |   4.00 |
          |  49.47 |   2.86 |
 ---------+--------+--------+
        1 |  28690 |  39703 |  68393
          |  29.38 |  40.66 |  70.04
          |  41.95 |  58.05 |
          |  50.53 |  97.14 |
 ---------+--------+--------+
 Total       56776    40874    97650
             58.14    41.86   100.00


TP= In both          = 39703
TN= In neither       = 28086
FP= In STAR only     = 1171
FN= In PacBio only   = 28690


FPR=FP/(FP+TN)=1171/(1171+28086)=4.00%
FNR=FN/(FN+TP)=28690/(28690+39703)=41.95%
*/

