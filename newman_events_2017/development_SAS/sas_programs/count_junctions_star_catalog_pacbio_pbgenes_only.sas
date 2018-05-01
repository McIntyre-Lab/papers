
/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* I am going to look at the false positive/negative rate for junction detection for STAR and Event Analysis
   by comparing these to the set of PacBio junctions

   False neg: detected with PB, but not with STAR/EA
   False pos: detected with STAR/EA, but not in PB

   Start with only "annotated junctions", "unannotated/novel junctions", then all. (skip border junctions)

   I am going to use a depth/average depth cut off of 5 reads, look at junctions detected in at least one rep, and then in both reps

   Note: FP rate will likely be very high for both STAR and EA. This is related to the level of depth of the PB experiment */

/* First: count junctions detected in each *
proc freq data=event.catalog_pacbio_star_junc_pb noprint;
   tables flag_in_pacbio*flag_in_catalog*flag_in_star / out=junc_by_set;
proc print data=junc_by_set;
run;

/*
       flag_in_    flag_in_    flag_in_
Obs     pacbio      catalog      star       COUNT    PERCENT

 1         0           0           1        19175     1.7705
 2         0           1           0       968402    89.4158
 3         0           1           1        22719     2.0977
 4         1           0           0         2755     0.2544
 5         1           0           1          936     0.0864
 6         1           1           0         2515     0.2322
 7         1           1           1        66530     6.1429

*/

proc freq data=event.catalog_pacbio_star_junc_pb;
   tables flag_in_pacbio*flag_in_star 
          flag_in_pacbio*flag_in_catalog  
          flag_in_catalog*flag_in_star ;
run;


/*
67466/109360 STAR junctions in PacBio (62%)
69045/1060166 catalog junctions in PacBio (7%)	-> expected: most of these are "unannotated"
89249/109360 STAR junctions in catalog (82%)
89249/1060166  catalog junctions in STAR (8%) 	-> expected: most of these are "unannotated"
*/

/* restrict to annotated junctions */

proc freq data=event.catalog_pacbio_star_junc_pb;
   where flag_junction_annotated=1;
   tables flag_in_pacbio*flag_in_star 
          flag_in_pacbio*flag_in_catalog  
          flag_in_catalog*flag_in_star ;
run;

/*
66117/82199 STAR junctions in PacBio (80%)
68383/97650 catalog junctions in PacBio (70%)
82199/82199 STAR junctions in catalog (100%)
82199/97650 catalog junctions in STAR (84%)
*/

/* Need to add a flag: is junction detected in BOTH samples (STAR/Event only) */

data flag_on_both;
   set event.catalog_pacbio_star_junc_pb;
   if flag_apn_NSC1_ge5=1 and flag_apn_NSC2_ge5=1 then flag_both_dtct_catalog=1;
   else flag_both_dtct_catalog=0;
   if flag_STAR_depth_NSC1_ge5=1 and flag_STAR_depth_NSC2_ge5=1 then flag_both_dtct_STAR=1;
   else flag_both_dtct_STAR=0;
run;

proc freq data=flag_on_both;
   where flag_in_star=1; *only count junctions in STAR;
   tables flag_both_dtct_STAR*flag_in_pacbio;
run;

/*
 flag_both_dtct_STAR
           flag_in_pacbio

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  38216 |  15236 |  53452
          |  34.95 |  13.93 |  48.88
          |  71.50 |  28.50 |
          |  91.22 |  22.58 |
 ---------+--------+--------+
        1 |   3678 |  52230 |  55908
          |   3.36 |  47.76 |  51.12
          |   6.58 |  93.42 |
          |   8.78 |  77.42 |
 ---------+--------+--------+
 Total       41894    67466   109360
             38.31    61.69   100.00

TP= in PB and dtct      = 52230
FP= not in PB and dtct  = 3678
TN= not in PB, not dtct = 38216
FN= in PB, not dtct     = 15236

FP rate = FP/(FP+TN) = 3678/(3678+38216) = 8.78%
FN rate = FN/(FN+TP) = 15236/(15236+52230) = 22.58%

*/

proc freq data=flag_on_both;
   where flag_intron_retention=0 and flag_in_catalog=1; *only count junctions that are in catalog;
   tables flag_both_dtct_catalog*flag_in_pacbio;
run;

/*


 flag_both_dtct_catalog
           flag_in_pacbio

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 905849 |  29311 | 935160
          |  92.80 |   3.00 |  95.80
          |  96.87 |   3.13 |
          |  99.86 |  42.45 |
 ---------+--------+--------+
        1 |   1233 |  39734 |  40967
          |   0.13 |   4.07 |   4.20
          |   3.01 |  96.99 |
          |   0.14 |  57.55 |
 ---------+--------+--------+
 Total      907082    69045   976127
             92.93     7.07   100.00

TP= in PB and dtct      = 39734
FP= not in PB and dtct  = 1233
TN= not in PB, not dtct = 905849
FN= in PB, not dtct     = 29311

FP rate = FP/(FP+TN) = 1233/(1233+905849) = 0.14 %
FN rate = FN/(FN+TP) = 29311/(29311+39734) = 42.45 %


FP rate here seems REALLY low ... BUT this is includes all the unannotated but plausible.
I should look at junctions that are annotated to known transcripts...
*/


/***** ANNOTATED JUNCTIONS ONLY *****/

proc freq data=flag_on_both;
   where flag_in_star=1 and flag_junction_annotated=1; *only count annotated junctions in STAR;
   tables flag_both_dtct_STAR*flag_in_pacbio;
run;

/*
 flag_both_dtct_STAR
           flag_in_pacbio

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  12889 |  14163 |  27052
          |  15.68 |  17.23 |  32.91
          |  47.65 |  52.35 |
          |  80.15 |  21.42 |
 ---------+--------+--------+
        1 |   3193 |  51954 |  55147
          |   3.88 |  63.21 |  67.09
          |   5.79 |  94.21 |
          |  19.85 |  78.58 |
 ---------+--------+--------+
 Total       16082    66117    82199
             19.56    80.44   100.00


TP= in PB and dtct      = 51954
FP= not in PB and dtct  = 3193
TN= not in PB, not dtct = 12889
FN= in PB, not dtct     = 14163

FP rate = FP/(FP+TN) = 3193/(3193+12889) = 19.85%
FN rate = FN/(FN+TP) = 14163/(14163+51954) = 21.42%
*/

proc freq data=flag_on_both;
   where flag_junction_annotated=1 and flag_in_catalog=1; *only count junctions that are in catalog;
   tables flag_both_dtct_catalog*flag_in_pacbio;
run;

/*
flag_both_dtct_catalog
          flag_in_pacbio

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  28086 |  28690 |  56776
         |  28.76 |  29.38 |  58.14
         |  49.47 |  50.53 |
         |  96.00 |  41.95 |
---------+--------+--------+
       1 |   1171 |  39703 |  40874
         |   1.20 |  40.66 |  41.86
         |   2.86 |  97.14 |
         |   4.00 |  58.05 |
---------+--------+--------+
Total       29257    68393    97650
            29.96    70.04   100.00

TP= in PB and dtct      = 39703
FP= not in PB and dtct  = 1171
TN= not in PB, not dtct = 28086
FN= in PB, not dtct     = 28690

FP rate = FP/(FP+TN) = 1171/(1171+28086) = 4.00 %
FN rate = FN/(FN+TP) = 28690/(28690+39703) = 41.95 %
*/


/***** UNANNOTATED/NOVEL EVENT ONLY *****/

proc freq data=flag_on_both;
   where flag_in_star=1 and flag_junction_annotated ne 1; *only count annotated junctions in STAR;
   tables flag_both_dtct_STAR*flag_in_pacbio;
run;

/*
flag_both_dtct_STAR
          flag_in_pacbio

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  25327 |   1073 |  26400
         |  93.25 |   3.95 |  97.20
         |  95.94 |   4.06 |
         |  98.12 |  79.54 |
---------+--------+--------+
       1 |    485 |    276 |    761
         |   1.79 |   1.02 |   2.80
         |  63.73 |  36.27 |
         |   1.88 |  20.46 |
---------+--------+--------+
Total       25812     1349    27161
            95.03     4.97   100.00



TP= in PB and dtct      = 276
FP= not in PB and dtct  = 485
TN= not in PB, not dtct = 25327
FN= in PB, not dtct     = 1073

FP rate = FP/(FP+TN) = 485/(485+25327) = 1.88%
FN rate = FN/(FN+TP) = 1073/(1073+276) = 79.54%
*/

proc freq data=flag_on_both;
   where flag_junction_annotated=0 and flag_intron_retention=0 and flag_in_catalog=1; *only count junctions that are in catalog;
   tables flag_both_dtct_catalog*flag_in_pacbio;
run;

/*

flag_both_dtct_catalog
          flag_in_pacbio

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 877763 |    621 | 878384
         |  99.92 |   0.07 |  99.99
         |  99.93 |   0.07 |
         |  99.99 |  95.25 |
---------+--------+--------+
       1 |     62 |     31 |     93
         |   0.01 |   0.00 |   0.01
         |  66.67 |  33.33 |
         |   0.01 |   4.75 |
---------+--------+--------+
Total      877825      652   878477
            99.93     0.07   100.00


TP= in PB and dtct      = 31
FP= not in PB and dtct  = 62
TN= not in PB, not dtct = 877763
FN= in PB, not dtct     = 621

FP rate = FP/(FP+TN) = 62/(62+877763) = 0.01 %
FN rate = FN/(FN+TP) = 621/(621+31) = 95.25 %

*/

