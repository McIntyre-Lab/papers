
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
   set event.catalog_pacbio_star_junc_pb;
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
        0 |  27744 |   4954 |  32698
          |  25.37 |   4.53 |  29.90
          |  84.85 |  15.15 |
          |  66.22 |   7.34 |
 ---------+--------+--------+
        1 |  14150 |  62512 |  76662
          |  12.94 |  57.16 |  70.10
          |  18.46 |  81.54 |
          |  33.78 |  92.66 |
 ---------+--------+--------+
 Total       41894    67466   109360
             38.31    61.69   100.00

TP= in PB and dtct      = 62512
FP= not in PB and dtct  = 14150
TN= not in PB, not dtct = 27744
FN= in PB, not dtct     = 4954

FP rate = FP/(FP+TN) = 14150/(14150+27744) = 33.78%
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
        0 | 896828 |   8743 | 905571
          |  91.88 |   0.90 |  92.77
          |  99.03 |   0.97 |
          |  98.87 |  12.66 |
 ---------+--------+--------+
        1 |  10254 |  60302 |  70556
          |   1.05 |   6.18 |   7.23
          |  14.53 |  85.47 |
          |   1.13 |  87.34 |
 ---------+--------+--------+
 Total      907082    69045   976127
             92.93     7.07   100.00

TP= in PB and dtct      = 60302
FP= not in PB and dtct  = 10254
TN= not in PB, not dtct = 896828
FN= in PB, not dtct     = 8743

FP rate = FP/(FP+TN) = 10254/(10254+896828) = 1.13 %
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
         0 |   6083 |   4367 |  10450
           |   7.40 |   5.31 |  12.71
           |  58.21 |  41.79 |
           |  37.82 |   6.60 |
  ---------+--------+--------+
         1 |   9999 |  61750 |  71749
           |  12.16 |  75.12 |  87.29
           |  13.94 |  86.06 |
           |  62.18 |  93.40 |
  ---------+--------+--------+
  Total       16082    66117    82199
              19.56    80.44   100.00


TP= in PB and dtct      = 61750
FP= not in PB and dtct  = 9999
TN= not in PB, not dtct = 6083
FN= in PB, not dtct     = 4367

FP rate = FP/(FP+TN) = 9999/(9999+6083) = 62.18%
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
         0 |  20654 |   8290 |  28944
           |  21.15 |   8.49 |  29.64
           |  71.36 |  28.64 |
           |  70.60 |  12.12 |
  ---------+--------+--------+
         1 |   8603 |  60103 |  68706
           |   8.81 |  61.55 |  70.36
           |  12.52 |  87.48 |
           |  29.40 |  87.88 |
  ---------+--------+--------+
  Total       29257    68393    97650
              29.96    70.04   100.00


TP= in PB and dtct      = 60103
FP= not in PB and dtct  = 8603
TN= not in PB, not dtct = 20654
FN= in PB, not dtct     = 8290

FP rate = FP/(FP+TN) = 8603/(8603+20654) = 29.40 %
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
         0 |  21661 |    587 |  22248
           |  79.75 |   2.16 |  81.91
           |  97.36 |   2.64 |
           |  83.92 |  43.51 |
  ---------+--------+--------+
         1 |   4151 |    762 |   4913
           |  15.28 |   2.81 |  18.09
           |  84.49 |  15.51 |
           |  16.08 |  56.49 |
  ---------+--------+--------+
  Total       25812     1349    27161
              95.03     4.97   100.00


TP= in PB and dtct      = 762
FP= not in PB and dtct  = 4151
TN= not in PB, not dtct = 21661
FN= in PB, not dtct     = 587

FP rate = FP/(FP+TN) = 4151/(4151+21661) = 16.08%
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
        0 | 876174 |    453 | 876627
          |  99.74 |   0.05 |  99.79
          |  99.95 |   0.05 |
          |  99.81 |  69.48 |
 ---------+--------+--------+
        1 |   1651 |    199 |   1850
          |   0.19 |   0.02 |   0.21
          |  89.24 |  10.76 |
          |   0.19 |  30.52 |
 ---------+--------+--------+
 Total      877825      652   878477
             99.93     0.07   100.00



TP= in PB and dtct      = 199
FP= not in PB and dtct  = 1651
TN= not in PB, not dtct = 876174
FN= in PB, not dtct     = 453

FP rate = FP/(FP+TN) = 1651/(1651+876174) = 0.19 %
FN rate = FN/(FN+TP) = 453/(453+199) = 69.48 %

*/

