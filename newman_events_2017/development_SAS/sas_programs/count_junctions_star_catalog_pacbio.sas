
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
proc freq data=event.catalog_pacbio_star_junctions noprint;
   tables flag_in_pacbio*flag_in_catalog*flag_in_star / out=junc_by_set;
proc print data=junc_by_set;
run;

/*

        flag_in_    flag_in_    flag_in_
 Obs     pacbio      catalog      star        COUNT    PERCENT

  1         0           0           1         49560     1.6240
  2         0           1           0       2848450    93.3396
  3         0           1           1         80959     2.6529
  4         1           0           0          2755     0.0903
  5         1           0           1           936     0.0307
  6         1           1           0          2515     0.0824
  7         1           1           1         66530     2.1801

*/

proc freq data=event.catalog_pacbio_star_junctions;
   tables flag_in_pacbio*flag_in_star 
          flag_in_pacbio*flag_in_catalog  
          flag_in_catalog*flag_in_star ;
run;


/*
67466/197985 STAR junctions in PacBio (34%)
69045/2998454 catalog junctions in PacBio (2%)	-> expected: monst of these are "unannotated"
147489/197985 STAR junctions in catalog (75%)
147489/2998454  catalog junctions in STAR (5%) 	-> expected: monst of these are "unannotated"
*/

/* restrict to annotated junctions */

proc freq data=event.catalog_pacbio_star_junctions;
   where flag_junction_annotated=1;
   tables flag_in_pacbio*flag_in_star 
          flag_in_pacbio*flag_in_catalog  
          flag_in_catalog*flag_in_star ;
run;

/*
66117/137925 STAR junctions in PacBio (48%)
68383/290055 catalog junctions in PacBio (24%)
137925/137925 STAR junctions in catalog (100%)
137925/290055 catalog junctions in STAR (48%)
*/

/* Need to add a flag: is junction detected in BOTH samples (STAR/Event only) */

data flag_on_both;
   set event.catalog_pacbio_star_junctions;
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
       0 | 114742 |  15236 | 129978
         |  57.95 |   7.70 |  65.65
         |  88.28 |  11.72 |
         |  87.91 |  22.58 |
---------+--------+--------+
       1 |  15777 |  52230 |  68007
         |   7.97 |  26.38 |  34.35
         |  23.20 |  76.80 |
         |  12.09 |  77.42 |
---------+--------+--------+
Total      130519    67466   197985
            65.92    34.08   100.00

TP= in PB and dtct      = 52230
FP= not in PB and dtct  = 15777
TN= not in PB, not dtct = 114742
FN= in PB, not dtct     = 15236

FP rate = FP/(FP+TN) = 15777/(15777+114742) = 12.09%
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
       0 |2674782 |  29311 |2704093
         |  97.25 |   1.07 |  98.31
         |  98.92 |   1.08 |
         |  99.75 |  42.45 |
---------+--------+--------+
       1 |   6644 |  39734 |  46378
         |   0.24 |   1.44 |   1.69
         |  14.33 |  85.67 |
         |   0.25 |  57.55 |
---------+--------+--------+
Total     2681426    69045  2750471
            97.49     2.51   100.00

TP= in PB and dtct      = 39734
FP= not in PB and dtct  = 6644
TN= not in PB, not dtct = 2674782
FN= in PB, not dtct     = 29311

FP rate = FP/(FP+TN) = 6644/(6644+2674782) = 0.25 %
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
        0 |  56932 |  14163 |  71095
          |  41.28 |  10.27 |  51.55
          |  80.08 |  19.92 |
          |  79.28 |  21.42 |
 ---------+--------+--------+
        1 |  14876 |  51954 |  66830
          |  10.79 |  37.67 |  48.45
          |  22.26 |  77.74 |
          |  20.72 |  78.58 |
 ---------+--------+--------+
 Total       71808    66117   137925
             52.06    47.94   100.00

TP= in PB and dtct      = 51954
FP= not in PB and dtct  = 14876
TN= not in PB, not dtct = 56932
FN= in PB, not dtct     = 14163

FP rate = FP/(FP+TN) = 14876/(14876+56932) = 20.72%
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
         0 | 215134 |  28690 | 243824
           |  74.17 |   9.89 |  84.06
           |  88.23 |  11.77 |
           |  97.05 |  41.95 |
  ---------+--------+--------+
         1 |   6528 |  39703 |  46231
           |   2.25 |  13.69 |  15.94
           |  14.12 |  85.88 |
           |   2.95 |  58.05 |
  ---------+--------+--------+
  Total      221662    68393   290055
              76.42    23.58   100.00


TP= in PB and dtct      = 39703
FP= not in PB and dtct  = 6528
TN= not in PB, not dtct = 215134
FN= in PB, not dtct     = 28690

FP rate = FP/(FP+TN) = 6528/(6528+215134) = 2.95 %
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
        0 |  57810 |   1073 |  58883
          |  96.25 |   1.79 |  98.04
          |  98.18 |   1.82 |
          |  98.47 |  79.54 |
 ---------+--------+--------+
        1 |    901 |    276 |   1177
          |   1.50 |   0.46 |   1.96
          |  76.55 |  23.45 |
          |   1.53 |  20.46 |
 ---------+--------+--------+
 Total       58711     1349    60060
             97.75     2.25   100.00


TP= in PB and dtct      = 276
FP= not in PB and dtct  = 901
TN= not in PB, not dtct = 57810
FN= in PB, not dtct     = 1073

FP rate = FP/(FP+TN) = 901/(901+57810) = 1.53%
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
        0 |2706071 |    621 |2706692
          |  99.91 |   0.02 |  99.94
          |  99.98 |   0.02 |
          |  99.94 |  95.25 |
 ---------+--------+--------+
        1 |   1676 |     31 |   1707
          |   0.06 |   0.00 |   0.06
          |  98.18 |   1.82 |
          |   0.06 |   4.75 |
 ---------+--------+--------+
 Total     2707747      652  2708399
             99.98     0.02   100.00



TP= in PB and dtct      = 31
FP= not in PB and dtct  = 1676
TN= not in PB, not dtct = 2706071
FN= in PB, not dtct     = 621

FP rate = FP/(FP+TN) = 1676/(1676+2706071) = 0.06 %
FN rate = FN/(FN+TP) = 621/(621+31) = 95.25 %

This doesn't seem right... I should limit this to ONLY junctions from genes that have PB transcripts
as PB is our "gold standard" and if a gene was not sequenced then we can't say much about it
*/

