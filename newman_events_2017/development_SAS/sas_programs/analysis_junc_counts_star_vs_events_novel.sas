
proc sort data=pb_nic_junc;
   by junc_index;
proc sort data=pb_nnc_junc;
   by junc_index;
proc sort data=star_novel;
   by junc_index;
proc sort data=star_unannot;
   by junc_index;
proc sort data=cat_unannot;
   by junc_index;
proc sort data=detection2;
   by junc_index;
run;

data star2event_novel;
  merge star_novel (in=in1) star_unannot (in=in2) cat_unannot (in=in3) pb_nic_junc (in=in4)
        pb_nnc_junc (in=in5) detection2;
  by junc_index;
  if in4 or in5 then flag_pb_junc=1; else flag_pb_junc=0;
  if in1 or in2 then flag_star_junc=1; else flag_star_junc=0;
  if in3 then flag_cat_junc=1; else flag_cat_junc=0;
  if in1 then flag_star_nnc=1; else flag_star_nnc=0;
  if in1 or in2 or in3 then output;
run;

* Check overlap;

*no NNC;
proc freq data=star2event_novel;
   where flag_star_nnc=0;
   tables flag_cat_junc*flag_star_junc;
run;

/*

  flag_cat_junc     flag_star_junc

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         1 |2698835 |   9564 |2708399
           |  99.65 |   0.35 | 100.00
           |  99.65 |   0.35 |
           | 100.00 | 100.00 |
  ---------+--------+--------+
  Total     2698835     9564  2708399
              99.65     0.35   100.00

9564/9564 (100%) unannotated STAR junctions are also unannotated catalog junctions (expected)
9564/2698835 (0.35%) unannotated catalog junctions are also in STAR

*/

*include NNC;
proc freq data=star2event_novel;
   tables flag_cat_junc*flag_star_junc;
run;


/*
 flag_cat_junc     flag_star_junc

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |  50496 |  50496
          |   0.00 |   1.83 |   1.83
          |   0.00 | 100.00 |
          |   0.00 |  84.08 |
 ---------+--------+--------+
        1 |2698835 |   9564 |2708399
          |  97.82 |   0.35 |  98.17
          |  99.65 |   0.35 |
          | 100.00 |  15.92 |
 ---------+--------+--------+
 Total     2698835    60060  2758895
             97.82     2.18   100.00

As above.
*/

/* NIC only */

* All;
proc freq data=star2event_novel;
   where flag_star_nnc=0;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*
flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |2689797 |  16353 |2706150
         |  99.31 |   0.60 |  99.92
         |  99.40 |   0.60 |
         |  99.96 |  92.79 |
---------+--------+--------+
       1 |    978 |   1271 |   2249
         |   0.04 |   0.05 |   0.08
         |  43.49 |  56.51 |
         |   0.04 |   7.21 |
---------+--------+--------+
Total     2690775    17624  2708399
            99.35     0.65   100.00

TP = In both    = 1271
TN = In neither = 2689797
FP = EA only    = 16353
FN = STAR only  = 978

FPR=FP/(FP+TN)=16353/(16353+2689797)= 0.60%
FNR=FN/(FN+TP)=978/(978+1271)= 43.49%


 flag_NPC_both_depth5
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |2706270 |   1609 |2707879
          |  99.92 |   0.06 |  99.98
          |  99.94 |   0.06 |
          |  99.98 |  94.26 |
 ---------+--------+--------+
        1 |    422 |     98 |    520
          |   0.02 |   0.00 |   0.02
          |  81.15 |  18.85 |
          |   0.02 |   5.74 |
 ---------+--------+--------+
 Total     2706692     1707  2708399
             99.94     0.06   100.00


TP = In both    = 98
TN = In neither = 2706270
FP = EA only    = 1609
FN = STAR only  = 422

FPR=FP/(FP+TN)=1609/(1609+2706270)=0.06%
FNR=FN/(FN+TP)=1079/(422+98)=81.15%


*/

* PB only;
proc freq data=star2event_novel;
   where flag_pb_junc=1 and flag_star_nnc=0;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*
 flag_NPC_both_depth0
           flag_NPC_both_apn0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    352 |     24 |    376
          |  53.99 |   3.68 |  57.67
          |  93.62 |   6.38 |
          |  77.70 |  12.06 |
 ---------+--------+--------+
        1 |    101 |    175 |    276
          |  15.49 |  26.84 |  42.33
          |  36.59 |  63.41 |
          |  22.30 |  87.94 |
 ---------+--------+--------+
 Total         453      199      652
             69.48    30.52   100.00

TP = In both    = 175
TN = In neither = 352
FP = EA only    = 24
FN = STAR only  = 101

FPR=FP/(FP+TN)=24/(24+352)= 6.38%
FNR=FN/(FN+TP)=101/(101+175)= 36.59%

flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |    522 |      1 |    523
         |  80.06 |   0.15 |  80.21
         |  99.81 |   0.19 |
         |  84.06 |   3.23 |
---------+--------+--------+
       1 |     99 |     30 |    129
         |  15.18 |   4.60 |  19.79
         |  76.74 |  23.26 |
         |  15.94 |  96.77 |
---------+--------+--------+
Total         621       31      652
            95.25     4.75   100.00


TP = In both    = 30
TN = In neither = 522
FP = EA only    = 1
FN = STAR only  = 99

FPR=FP/(FP+TN)=1/(1+522)= 0.19%
FNR=FN/(FN+TP)=99/(99+30)= 76.74%

*/

* Gene has PB;
proc freq data=star2event_novel;
   where flag_gene_has_pb=1 and flag_star_nnc=0;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;

/*

flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 950007 |  10789 | 960796
         |  98.70 |   1.12 |  99.82
         |  98.88 |   1.12 |
         |  99.92 |  91.72 |
---------+--------+--------+
       1 |    746 |    974 |   1720
         |   0.08 |   0.10 |   0.18
         |  43.37 |  56.63 |
         |   0.08 |   8.28 |
---------+--------+--------+
Total      950753    11763   962516
            98.78     1.22   100.00

TP = In both    = 974
TN = In neither = 950007
FP = EA only    = 10789
FN = STAR only  = 746

FPR=FP/(FP+TN)=10789/(10789+950007)= 1.12%
FNR=FN/(FN+TP)=746/(746+974)= 43.37%

 flag_NPC_both_depth5
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 961329 |    787 | 962116
          |  99.88 |   0.08 |  99.96
          |  99.92 |   0.08 |
          |  99.96 |  92.59 |
 ---------+--------+--------+
        1 |    337 |     63 |    400
          |   0.04 |   0.01 |   0.04
          |  84.25 |  15.75 |
          |   0.04 |   7.41 |
 ---------+--------+--------+
 Total      961666      850   962516
             99.91     0.09   100.00


TP = In both    = 63
TN = In neither = 961329
FP = EA only    = 787
FN = STAR only  = 337

FPR=FP/(FP+TN)=787/(787+961329)= 0.08%
FNR=FN/(FN+TP)=337/(337+63)= 84.25%
*/



/* Unannot EA vs NIC/NNC STAR */

* All;
proc freq data=star2event_novel;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*
   flag_NPC_both_depth0
             flag_NPC_both_apn0

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |2734891 |  16353 |2751244
            |  99.13 |   0.59 |  99.72
            |  99.41 |   0.59 |
            |  99.77 |  92.79 |
   ---------+--------+--------+
          1 |   6380 |   1271 |   7651
            |   0.23 |   0.05 |   0.28
            |  83.39 |  16.61 |
            |   0.23 |   7.21 |
   ---------+--------+--------+
   Total     2741271    17624  2758895
               99.36     0.64   100.00

TP = In both    = 1271
TN = In neither = 2734891
FP = EA only    = 16353
FN = STAR only  = 6380

FPR=FP/(FP+TN)=16353/(16353+2734891)= 0.59%
FNR=FN/(FN+TP)=6380/(6380+1271)= 83.39%


flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |2756109 |   1609 |2757718
         |  99.90 |   0.06 |  99.96
         |  99.94 |   0.06 |
         |  99.96 |  94.26 |
---------+--------+--------+
       1 |   1079 |     98 |   1177
         |   0.04 |   0.00 |   0.04
         |  91.67 |   8.33 |
         |   0.04 |   5.74 |
---------+--------+--------+
Total     2757188     1707  2758895
            99.94     0.06   100.00


TP = In both    = 98
TN = In neither = 2756109
FP = EA only    = 1609
FN = STAR only  = 1079

FPR=FP/(FP+TN)=1609/(1609+2756109)=0.06%
FNR=FN/(FN+TP)=1079/(1079+98)=91.67%


*/

* PB only;
proc freq data=star2event_novel;
   where flag_pb_junc=1;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;


/*

 flag_NPC_both_depth0
           flag_NPC_both_apn0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    802 |     24 |    826
          |  50.50 |   1.51 |  52.02
          |  97.09 |   2.91 |
          |  57.74 |  12.06 |
 ---------+--------+--------+
        1 |    587 |    175 |    762
          |  36.96 |  11.02 |  47.98
          |  77.03 |  22.97 |
          |  42.26 |  87.94 |
 ---------+--------+--------+
 Total        1389      199     1588
             87.47    12.53   100.00

TP = In both    = 175
TN = In neither = 802
FP = EA only    = 24
FN = STAR only  = 587

FPR=FP/(FP+TN)=24/(24+802)= 2.91%
FNR=FN/(FN+TP)=587/(587+175)= 77.03%


 flag_NPC_both_depth5
           flag_NPC_both_apn5

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   1311 |      1 |   1312
          |  82.56 |   0.06 |  82.62
          |  99.92 |   0.08 |
          |  84.20 |   3.23 |
 ---------+--------+--------+
        1 |    246 |     30 |    276
          |  15.49 |   1.89 |  17.38
          |  89.13 |  10.87 |
          |  15.80 |  96.77 |
 ---------+--------+--------+
 Total        1557       31     1588
             98.05     1.95   100.00


TP = In both    = 30
TN = In neither = 1311
FP = EA only    = 1
FN = STAR only  = 246

FPR=FP/(FP+TN)=1/(1+1311)= 0.08%
FNR=FN/(FN+TP)=246/(246+30)= 89.13%

*/

* Gene has PB;
proc freq data=star2event_novel;
   where flag_gene_has_pb=1;
   tables flag_NPC_both_depth0*flag_NPC_both_apn0
          flag_NPC_both_depth5*flag_NPC_both_apn5;
run;

/*

flag_NPC_both_depth0
          flag_NPC_both_apn0

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 966925 |  10789 | 977714
         |  98.40 |   1.10 |  99.50
         |  98.90 |   1.10 |
         |  99.59 |  91.72 |
---------+--------+--------+
       1 |   3939 |    974 |   4913
         |   0.40 |   0.10 |   0.50
         |  80.18 |  19.82 |
         |   0.41 |   8.28 |
---------+--------+--------+
Total      970864    11763   982627
            98.80     1.20   100.00

TP = In both    = 974
TN = In neither = 966925
FP = EA only    = 10789
FN = STAR only  = 3939

FPR=FP/(FP+TN)=10789/(10789+966925)= 1.10%
FNR=FN/(FN+TP)=3939/(3939+974)= 80.18%


flag_NPC_both_depth5
          flag_NPC_both_apn5

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 981079 |    787 | 981866
         |  99.84 |   0.08 |  99.92
         |  99.92 |   0.08 |
         |  99.93 |  92.59 |
---------+--------+--------+
       1 |    698 |     63 |    761
         |   0.07 |   0.01 |   0.08
         |  91.72 |   8.28 |
         |   0.07 |   7.41 |
---------+--------+--------+
Total      981777      850   982627
            99.91     0.09   100.00


TP = In both    = 63
TN = In neither = 981079
FP = EA only    = 787
FN = STAR only  = 698

FPR=FP/(FP+TN)=787/(787+981079)= 0.08%
FNR=FN/(FN+TP)=698/(698+63)= 91.72%
*/


