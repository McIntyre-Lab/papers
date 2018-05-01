
/* Compare EA to PacBio -- unannotated junctions */

proc sort data=pb_nic_junc;
   by junc_index;
proc sort data=cat_unannot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data pb2cat_nic;
  merge pb_nic_junc (in=in1) cat_unannot (in=in2) detection2 (in=in3);
  by junc_index;
  if in1 then flag_junc_in_pacbio=1; else flag_junc_in_pacbio=0;
  if in2 then flag_junc_in_catalog=1; else flag_junc_in_catalog=0;
  if in1 or in2;
run;


proc freq data=pb2cat_nic;
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
        0 |      0 |  52391 |  52391
          |   0.00 |  98.77 |  98.77
          |   0.00 | 100.00 |
          |   0.00 |  99.38 |
 ---------+--------+--------+
        1 |    325 |    327 |    652
          |   0.61 |   0.62 |   1.23
          |  49.85 |  50.15 |
          | 100.00 |   0.62 |
 ---------+--------+--------+
 Total         325    52718    53043
              0.61    99.39   100.00



50.15% of NIC PacBio junctions are "detectable" (at least one read aligns)
Though, these only make up 0.62% of all detectable unannotated junctions

*/

* all;
proc freq data=pb2cat_nic;
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
        0 |  34966 |  17425 |  52391
          |  65.92 |  32.85 |  98.77
          |  66.74 |  33.26 |
          |  98.72 |  98.87 |
 ---------+--------+--------+
        1 |    453 |    199 |    652
          |   0.85 |   0.38 |   1.23
          |  69.48 |  30.52 |
          |   1.28 |   1.13 |
 ---------+--------+--------+
 Total       35419    17624    53043
             66.77    33.23   100.00


TP= In both          = 199
TN= In neither       = 34966
FP= In STAR only     = 17425
FN= In PacBio only   = 453

FPR=FP/(FP+TN)= 17425/(17425+34966)= 33.26%
FNR=FN/(FN+TP)= 453/(453+199)= 69.48%

flag_junc_in_pacbio
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  50715 |   1676 |  52391
         |  95.61 |   3.16 |  98.77
         |  96.80 |   3.20 |
         |  98.79 |  98.18 |
---------+--------+--------+
       1 |    621 |     31 |    652
         |   1.17 |   0.06 |   1.23
         |  95.25 |   4.75 |
         |   1.21 |   1.82 |
---------+--------+--------+
Total       51336     1707    53043
            96.78     3.22   100.00




TP= In both          = 31
TN= In neither       = 50715
FP= In STAR only     = 1676
FN= In PacBio only   = 621

FPR=FP/(FP+TN)= 1676/(1676+2706071)= 3.20%
FNR=FN/(FN+TP)= 621/(621+31)= 95.25%
*/

* PB junc only;
proc freq data=pb2cat_nic;
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
        1 |    453 |    199 |    652
          |  69.48 |  30.52 | 100.00
          |  69.48 |  30.52 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total         453      199      652
             69.48    30.52   100.00


TP= In both          = 199Unannotated/”Novel in catalog” junctions: STAR vs EA
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 453


FPR=FP/(FP+TN)= n/a
FNR=FN/(FN+TP)= 453/(453+199)= 69.48%

 flag_junc_in_pacbio
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        1 |    621 |     31 |    652
          |  95.25 |   4.75 | 100.00
          |  95.25 |   4.75 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total         621       31      652
             95.25     4.75   100.00


TP= In both          = 31
TN= In neither       = 0
FP= In STAR only     = 0
FN= In PacBio only   = 621

FPR=FP/(FP+TN)=n/a
FNR=FN/(FN+TP)= 621/(621+31)= 95.25%
*/

* Gene has PB;
proc freq data=pb2cat_nic;
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
       0 |  22152 |  11564 |  33716
         |  64.46 |  33.65 |  98.10
         |  65.70 |  34.30 |
         |  98.00 |  98.31 |
---------+--------+--------+
       1 |    453 |    199 |    652
         |   1.32 |   0.58 |   1.90
         |  69.48 |  30.52 |
         |   2.00 |   1.69 |
---------+--------+--------+
Total       22605    11763    34368
            65.77    34.23   100.00


TP= In both          = 199
TN= In neither       = 22152
FP= In STAR only     = 11564
FN= In PacBio only   = 453


FPR=FP/(FP+TN)=11564/(11564+22152)=34.30%
FNR=FN/(FN+TP)=453/(453+199)=69.48%


flag_junc_in_pacbio
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  32897 |    819 |  33716
         |  95.72 |   2.38 |  98.10
         |  97.57 |   2.43 |
         |  98.15 |  96.35 |
---------+--------+--------+
       1 |    621 |     31 |    652
         |   1.81 |   0.09 |   1.90
         |  95.25 |   4.75 |
         |   1.85 |   3.65 |
---------+--------+--------+
Total       33518      850    34368
            97.53     2.47   100.00


TP= In both          = 31
TN= In neither       = 32897
FP= In STAR only     = 819
FN= In PacBio only   = 621


FPR=FP/(FP+TN)=819/(819+32897)=2.43%
FNR=FN/(FN+TP)=621/(621+31)=95.25%
*/

