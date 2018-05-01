
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
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |  63015 |  63015
          |   0.00 |  47.95 |  47.95
          |   0.00 | 100.00 |
          |   0.00 |  49.64 |
 ---------+--------+--------+
        1 |   4473 |  63920 |  68393
          |   3.40 |  48.64 |  52.05
          |   6.54 |  93.46 |
          | 100.00 |  50.36 |
 ---------+--------+--------+
 Total        4473   126935   131408
              3.40    96.60   100.00

93.46% of annotated PacBio junctions are in the junction catalog
Though, these make up 50.36% of all detectable annotated junctions

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
          0 |  27736 |  35279 |  63015
            |  21.11 |  26.85 |  47.95
            |  44.01 |  55.99 |
            |  76.99 |  36.99 |
   ---------+--------+--------+
          1 |   8290 |  60103 |  68393
            |   6.31 |  45.74 |  52.05
            |  12.12 |  87.88 |
            |  23.01 |  63.01 |
   ---------+--------+--------+
   Total       36026    95382   131408
               27.42    72.58   100.00




TP= In both          = 60103
TN= In neither       = 27736
FP= In STAR only     = 35279
FN= In PacBio only   = 8290

FPR=FP/(FP+TN)= 35279/(35279+27736)= 55.99%
FNR=FN/(FN+TP)= 8290/(8290+60103)= 12.12%


flag_junc_in_pacbio
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  56487 |   6528 |  63015
         |  42.99 |   4.97 |  47.95
         |  89.64 |  10.36 |
         |  66.32 |  14.12 |
---------+--------+--------+
       1 |  28690 |  39703 |  68393
         |  21.83 |  30.21 |  52.05
         |  41.95 |  58.05 |
         |  33.68 |  85.88 |
---------+--------+--------+
Total       85177    46231   131408
            64.82    35.18   100.00



TP= In both          = 39703
TN= In neither       = 56487
FP= In STAR only     = 6528
FN= In PacBio only   = 28690

FPR=FP/(FP+TN)= 6528/(6528+56487)= 10.36%
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
         0 |   6258 |   8603 |  14861
           |   7.52 |  10.33 |  17.85
           |  42.11 |  57.89 |
           |  43.02 |  12.52 |
  ---------+--------+--------+
         1 |   8290 |  60103 |  68393
           |   9.96 |  72.19 |  82.15
           |  12.12 |  87.88 |
           |  56.98 |  87.48 |
  ---------+--------+--------+
  Total       14548    68706    83254
              17.47    82.53   100.00


TP= In both          = 60103
TN= In neither       = 6258
FP= In STAR only     = 8603
FN= In PacBio only   = 8290


FPR=FP/(FP+TN)=8603/(8603+6258)=57.89%
FNR=FN/(FN+TP)=8290/(8290+60103)=12.12%

  flag_junc_in_pacbio
            flag_NPC_both_apn5

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  13690 |   1171 |  14861
           |  16.44 |   1.41 |  17.85
           |  92.12 |   7.88 |
           |  32.30 |   2.86 |
  ---------+--------+--------+
         1 |  28690 |  39703 |  68393
           |  34.46 |  47.69 |  82.15
           |  41.95 |  58.05 |
           |  67.70 |  97.14 |
  ---------+--------+--------+
  Total       42380    40874    83254
              50.90    49.10   100.00



TP= In both          = 39703
TN= In neither       = 13690
FP= In STAR only     = 1171
FN= In PacBio only   = 28690

FPR=FP/(FP+TN)=1171/(1171+13690)= 7.88%
FNR=FN/(FN+TP)=28690/(28690+39703)=41.95%
*/


