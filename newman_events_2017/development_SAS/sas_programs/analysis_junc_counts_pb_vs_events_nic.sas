
/* Compare STAR to PacBio -- unannotated junctions */

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
Col Pct  |       1|  Total
---------+--------+
       0 |2707747 |2707747
         |  99.98 |  99.98
         | 100.00 |
         |  99.98 |
---------+--------+
       1 |    652 |    652
         |   0.02 |   0.02
         | 100.00 |
         |   0.02 |
---------+--------+
Total     2708399  2708399
           100.00   100.00


100% of NIC PacBio junctions are in the junction catalog (expected)
Though, these only make up 0.02% of all unannotated junctions

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
       0 |2690322 |  17425 |2707747
         |  99.33 |   0.64 |  99.98
         |  99.36 |   0.64 |
         |  99.98 |  98.87 |
---------+--------+--------+
       1 |    453 |    199 |    652
         |   0.02 |   0.01 |   0.02
         |  69.48 |  30.52 |
         |   0.02 |   1.13 |
---------+--------+--------+
Total     2690775    17624  2708399
            99.35     0.65   100.00


TP= In both          = 199
TN= In neither       = 2690322
FP= In STAR only     = 17425
FN= In PacBio only   = 453

FPR=FP/(FP+TN)= 17425/(17425+2690322)= 0.64%
FNR=FN/(FN+TP)= 453/(453+199)= 30.52%

 flag_junc_in_pacbio
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |2706071 |   1676 |2707747
          |  99.91 |   0.06 |  99.98
          |  99.94 |   0.06 |
          |  99.98 |  98.18 |
 ---------+--------+--------+
        1 |    621 |     31 |    652
          |   0.02 |   0.00 |   0.02
          |  95.25 |   4.75 |
          |   0.02 |   1.82 |
 ---------+--------+--------+
 Total     2706692     1707  2708399
             99.94     0.06   100.00

TP= In both          = 31
TN= In neither       = 2706071
FP= In STAR only     = 1676
FN= In PacBio only   = 621

FPR=FP/(FP+TN)= 1676/(1676+2706071)= 0.06%
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


TP= In both          = 199
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
       0 | 950300 |  11564 | 961864
         |  98.73 |   1.20 |  99.93
         |  98.80 |   1.20 |
         |  99.95 |  98.31 |
---------+--------+--------+
       1 |    453 |    199 |    652
         |   0.05 |   0.02 |   0.07
         |  69.48 |  30.52 |
         |   0.05 |   1.69 |
---------+--------+--------+
Total      950753    11763   962516
            98.78     1.22   100.00

           The SAS System         09

TP= In both          = 199
TN= In neither       = 950300
FP= In STAR only     = 11564
FN= In PacBio only   = 453


FPR=FP/(FP+TN)=11564/(11564+950300)=1.20%
FNR=FN/(FN+TP)=453/(453+199)=69.48%

flag_junc_in_pacbio
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 961045 |    819 | 961864
         |  99.85 |   0.09 |  99.93
         |  99.91 |   0.09 |
         |  99.94 |  96.35 |
---------+--------+--------+
       1 |    621 |     31 |    652
         |   0.06 |   0.00 |   0.07
         |  95.25 |   4.75 |
         |   0.06 |   3.65 |
---------+--------+--------+
Total      961666      850   962516
            99.91     0.09   100.00



TP= In both          = 31
TN= In neither       = 961045
FP= In STAR only     = 819
FN= In PacBio only   = 621


FPR=FP/(FP+TN)=819/(819+961045)=0.09%
FNR=FN/(FN+TP)=621/(621+31)=95.25%
*/

