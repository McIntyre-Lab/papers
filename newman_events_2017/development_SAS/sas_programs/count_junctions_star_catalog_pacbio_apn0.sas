
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

/* Need to add a flag: is junction detected in BOTH samples (STAR/Event only) */

data flag_on_both;
   set event.catalog_pacbio_star_junctions;
   if flag_apn_NSC1_gt0=1 and flag_apn_NSC2_gt0=1 then flag_both_dtct_catalog=1;
   else flag_both_dtct_catalog=0;
   if flag_STAR_depth_NSC1_gt0=1 and flag_STAR_depth_NSC2_gt0=1 then flag_both_dtct_STAR=1;
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
       0 |  81699 |   4954 |  86653
         |  41.27 |   2.50 |  43.77
         |  94.28 |   5.72 |
         |  62.60 |   7.34 |
---------+--------+--------+
       1 |  48820 |  62512 | 111332
         |  24.66 |  31.57 |  56.23
         |  43.85 |  56.15 |
         |  37.40 |  92.66 |
---------+--------+--------+
Total      130519    67466   197985
            65.92    34.08   100.00

TP= in PB and dtct      = 62512
FP= not in PB and dtct  = 48820
TN= not in PB, not dtct = 81699
FN= in PB, not dtct     = 4954

FP rate = FP/(FP+TN) = 48820/(48820+81699) = 37.40%
FN rate = FN/(FN+TP) = 4954/(4954+62512) = 7.34%

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
        0 |2643588 |   8743 |2652331
          |  96.11 |   0.32 |  96.43
          |  99.67 |   0.33 |
          |  98.59 |  12.66 |
 ---------+--------+--------+
        1 |  37838 |  60302 |  98140
          |   1.38 |   2.19 |   3.57
          |  38.56 |  61.44 |
          |   1.41 |  87.34 |
 ---------+--------+--------+
 Total     2681426    69045  2750471
             97.49     2.51   100.00



TP= in PB and dtct      = 60302
FP= not in PB and dtct  = 37838
TN= not in PB, not dtct = 2643588
FN= in PB, not dtct     = 8743

FP rate = FP/(FP+TN) = 37838/(37838+2643588) = 1.41 %
FN rate = FN/(FN+TP) = 8743/(8743+60302) = 12.66 %

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
        0 |  29877 |   4367 |  34244
          |  21.66 |   3.17 |  24.83
          |  87.25 |  12.75 |
          |  41.61 |   6.60 |
 ---------+--------+--------+
        1 |  41931 |  61750 | 103681
          |  30.40 |  44.77 |  75.17
          |  40.44 |  59.56 |
          |  58.39 |  93.40 |
 ---------+--------+--------+
 Total       71808    66117   137925
             52.06    47.94   100.00


TP= in PB and dtct      = 61750
FP= not in PB and dtct  = 41931
TN= not in PB, not dtct = 29877
FN= in PB, not dtct     = 4367

FP rate = FP/(FP+TN) = 41931/(41931+29877) = 58.39%
FN rate = FN/(FN+TP) = 4367/(4367+61750) = 6.60%
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
       0 | 186383 |   8290 | 194673
         |  64.26 |   2.86 |  67.12
         |  95.74 |   4.26 |
         |  84.08 |  12.12 |
---------+--------+--------+
       1 |  35279 |  60103 |  95382
         |  12.16 |  20.72 |  32.88
         |  36.99 |  63.01 |
         |  15.92 |  87.88 |
---------+--------+--------+
Total      221662    68393   290055
            76.42    23.58   100.00 

TP= in PB and dtct      = 60103
FP= not in PB and dtct  = 35279
TN= not in PB, not dtct = 186383
FN= in PB, not dtct     = 8290

FP rate = FP/(FP+TN) = 35279/(35279+186383) = 15.92 %
FN rate = FN/(FN+TP) = 8290/(8290+60103) = 12.12 %
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
        0 |  51822 |    587 |  52409
          |  86.28 |   0.98 |  87.26
          |  98.88 |   1.12 |
          |  88.27 |  43.51 |
 ---------+--------+--------+
        1 |   6889 |    762 |   7651
          |  11.47 |   1.27 |  12.74
          |  90.04 |   9.96 |
          |  11.73 |  56.49 |
 ---------+--------+--------+
 Total       58711     1349    60060
             97.75     2.25   100.00


TP= in PB and dtct      = 762
FP= not in PB and dtct  = 6889
TN= not in PB, not dtct = 51822
FN= in PB, not dtct     = 587

FP rate = FP/(FP+TN) = 6889/(6889+51822) = 11.73%
FN rate = FN/(FN+TP) = 587/(587+762) = 43.51%
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
        0 |2457205 |    453 |2457658
          |  99.87 |   0.02 |  99.89
          |  99.98 |   0.02 |
          |  99.90 |  69.48 |
 ---------+--------+--------+
        1 |   2559 |    199 |   2758
          |   0.10 |   0.01 |   0.11
          |  92.78 |   7.22 |
          |   0.10 |  30.52 |
 ---------+--------+--------+
 Total     2459764      652  2460416
             99.97     0.03   100.00

TP= in PB and dtct      = 199
FP= not in PB and dtct  = 2559
TN= not in PB, not dtct = 2457205
FN= in PB, not dtct     = 453

FP rate = FP/(FP+TN) = 2559/(2559+2457205) = 0.10 %
FN rate = FN/(FN+TP) = 453/(453+199) = 69.48 %

This doesn't seem right... I should limit this to ONLY junctions from genes that have PB transcripts
as PB is our "gold standard" and if a gene was not sequenced then we can't say much about it
*/

