ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Now I am going to set some flags to help me summarize the data.
   I want to flag if the within-group difference is larger than the between-group difference first.
   Then for the whole-group comparisons (ie, where I have merged reps),
   flag if diff is >0.2, >0.5 and >0.8, and flag if bayes factor is >5, >10, >20. */


/* Flag if within cell type difference is larger than between cell type difference */

data flag_within_gt_btwn_diffs;
  set event.miso_all_results_nsc_v_old_trim;
  /* Flag if within group delta-psi's are greater than between group delta-psi's */

  if NSC1_NSC2_diff ne . and NSC1_OLD1_diff ne . then do;
         if abs(NSC1_NSC2_diff) > abs(NSC1_OLD1_diff) then flag_diff_NSC1v2_gt_NSC1vOLD1=1;
         else flag_diff_NSC1v2_gt_NSC1vOLD1=0; end;

  if NSC1_NSC2_diff ne . and NSC1_OLD2_diff ne . then do;
         if abs(NSC1_NSC2_diff) > abs(NSC1_OLD2_diff) then flag_diff_NSC1v2_gt_NSC1vOLD2=1;
         else flag_diff_NSC1v2_gt_NSC1vOLD2=0; end;

  if NSC1_NSC2_diff ne . and NSC2_OLD1_diff ne . then do;
         if abs(NSC1_NSC2_diff) > abs(NSC2_OLD1_diff) then flag_diff_NSC1v2_gt_NSC2vOLD1=1;
         else flag_diff_NSC1v2_gt_NSC2vOLD1=0; end;

  if NSC1_NSC2_diff ne . and NSC2_OLD2_diff ne . then do;
         if abs(NSC1_NSC2_diff) > abs(NSC2_OLD2_diff) then flag_diff_NSC1v2_gt_NSC2vOLD2=1;
         else flag_diff_NSC1v2_gt_NSC2vOLD2=0; end;

  if OLD1_OLD2_diff ne . and NSC1_OLD1_diff ne . then do;
         if abs(OLD1_OLD2_diff) > abs(NSC1_OLD1_diff) then flag_diff_OLD1v2_gt_NSC1vOLD1=1;
         else flag_diff_OLD1v2_gt_NSC1vOLD1=0; end;

  if OLD1_OLD2_diff ne . and NSC1_OLD2_diff ne . then do;
         if abs(OLD1_OLD2_diff) > abs(NSC1_OLD2_diff) then flag_diff_OLD1v2_gt_NSC1vOLD2=1;
         else flag_diff_OLD1v2_gt_NSC1vOLD2=0; end;

  if OLD1_OLD2_diff ne . and NSC2_OLD1_diff ne . then do;
         if abs(OLD1_OLD2_diff) > abs(NSC2_OLD1_diff) then flag_diff_OLD1v2_gt_NSC2vOLD1=1;
         else flag_diff_OLD1v2_gt_NSC2vOLD1=0; end;

  if OLD1_OLD2_diff ne . and NSC2_OLD2_diff ne . then do;
         if abs(OLD1_OLD2_diff) > abs(NSC2_OLD2_diff) then flag_diff_OLD1v2_gt_NSC2vOLD2=1;
         else flag_diff_OLD1v2_gt_NSC2vOLD2=0; end;

   if flag_diff_NSC1v2_gt_NSC1vOLD1=1 or flag_diff_NSC1v2_gt_NSC1vOLD2=1
   or flag_diff_NSC1v2_gt_NSC2vOLD1=1 or flag_diff_NSC1v2_gt_NSC2vOLD2=1
   or flag_diff_OLD1v2_gt_NSC1vOLD1=1 or flag_diff_OLD1v2_gt_NSC1vOLD2=1
   or flag_diff_OLD1v2_gt_NSC2vOLD1=1 or flag_diff_OLD1v2_gt_NSC2vOLD2=1
   then flag_within_diff_gt_btwn_diff=1; else flag_within_diff_gt_btwn_diff=0;
run;

proc freq data=flag_within_gt_btwn_diffs;
   tables flag_diff_NSC1v2_gt_NSC1vOLD1 flag_diff_NSC1v2_gt_NSC1vOLD2
          flag_diff_NSC1v2_gt_NSC2vOLD1 flag_diff_NSC1v2_gt_NSC2vOLD2
          flag_diff_OLD1v2_gt_NSC1vOLD1 flag_diff_OLD1v2_gt_NSC1vOLD2
          flag_diff_OLD1v2_gt_NSC2vOLD1 flag_diff_OLD1v2_gt_NSC2vOLD2
          flag_within_diff_gt_btwn_diff;
run;

/* Summary:
2048 instances where within-NSC diff > NSC1 vs OLD1
2134 instances where within-NSC diff > NSC1 vs OLD2
1983 instances where within-NSC diff > NSC2 vs OLD1
2020 instances where within-NSC diff > NSC2 vs OLD2

1778 instances where within-OLD diff > NSC1 vs OLD1
2034 instances where within-OLD diff > NSC1 vs OLD2
1815 instances where within-OLD diff > NSC2 vs OLD1
1941 instances where within-OLD diff > NSC2 vs OLD2

In total 4337 instances where the within-group diff > between-group diff
*/

/* On the main comparison, I now want to bin the delta psi's and bayes factors */


data bin_diffs_and_bayes;
   set flag_within_gt_btwn_diffs;
   if abs(NSC_OLD_diff) = . then bin_NSCvOLD_diff=.;
   else if abs(NSC_OLD_diff) lt 0.2 then bin_NSCvOLD_diff=0; *abs diff < 0.2;
   else if abs(NSC_OLD_diff) lt 0.5 then bin_NSCvOLD_diff=1;  *abs diff < 0.5;
   else if abs(NSC_OLD_diff) lt 0.8 then bin_NSCvOLD_diff=2;  *abs diff < 0.8;
   else bin_NSCvOLD_diff=3;  *abs diff > 0.8;

   if NSC_OLD_bayes_factor =. then bin_NSCvOLD_bayes=.;
   else if NSC_OLD_bayes_factor lt 5 then bin_NSCvOLD_bayes=0;
   else if NSC_OLD_bayes_factor lt 10 then bin_NSCvOLD_bayes=1;
   else if NSC_OLD_bayes_factor lt 20 then bin_NSCvOLD_bayes=2;
   else bin_NSCvOLD_bayes=3;
run;

/* Count */

proc freq data=bin_diffs_and_bayes;
   tables bin_NSCvOLD_diff*bin_NSCvOLD_bayes;
run;

/*

bin_NSCvOLD_diff     bin_NSCvOLD_bayes

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|       2|       3|  Total
---------+--------+--------+--------+--------+
       0 |   5236 |     16 |     12 |     61 |   5325
         |  85.89 |   0.26 |   0.20 |   1.00 |  87.35
         |  98.33 |   0.30 |   0.23 |   1.15 |
         |  88.85 |  44.44 |  36.36 |  45.52 |
---------+--------+--------+--------+--------+
       1 |    653 |     17 |     19 |     57 |    746
         |  10.71 |   0.28 |   0.31 |   0.94 |  12.24
         |  87.53 |   2.28 |   2.55 |   7.64 |
         |  11.08 |  47.22 |  57.58 |  42.54 |
---------+--------+--------+--------+--------+
       2 |      4 |      3 |      2 |     15 |     24
         |   0.07 |   0.05 |   0.03 |   0.25 |   0.39
         |  16.67 |  12.50 |   8.33 |  62.50 |
         |   0.07 |   8.33 |   6.06 |  11.19 |
---------+--------+--------+--------+--------+
       3 |      0 |      0 |      0 |      1 |      1
         |   0.00 |   0.00 |   0.00 |   0.02 |   0.02
         |   0.00 |   0.00 |   0.00 | 100.00 |
         |   0.00 |   0.00 |   0.00 |   0.75 |
---------+--------+--------+--------+--------+
Total        5893       36       33      134     6096
            96.67     0.59     0.54     2.20   100.00

Frequency Missing = 242
*/

/* And for only events that are not different within groups */

proc freq data=bin_diffs_and_bayes;
   where flag_within_diff_gt_btwn_diff=0;
   tables bin_NSCvOLD_diff*bin_NSCvOLD_bayes;
run;

/*
    Table of bin_NSCvOLD_diff by bin_NSCvOLD_bayes

bin_NSCvOLD_diff     bin_NSCvOLD_bayes

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|       2|       3|  Total
---------+--------+--------+--------+--------+
       0 |   1309 |      7 |      7 |     42 |   1365
         |  74.42 |   0.40 |   0.40 |   2.39 |  77.60
         |  95.90 |   0.51 |   0.51 |   3.08 |
         |  80.90 |  36.84 |  31.82 |  42.00 |
---------+--------+--------+--------+--------+
       1 |    305 |     11 |     13 |     44 |    373
         |  17.34 |   0.63 |   0.74 |   2.50 |  21.21
         |  81.77 |   2.95 |   3.49 |  11.80 |
         |  18.85 |  57.89 |  59.09 |  44.00 |
---------+--------+--------+--------+--------+
       2 |      4 |      1 |      2 |     13 |     20
         |   0.23 |   0.06 |   0.11 |   0.74 |   1.14
         |  20.00 |   5.00 |  10.00 |  65.00 |
         |   0.25 |   5.26 |   9.09 |  13.00 |
---------+--------+--------+--------+--------+
       3 |      0 |      0 |      0 |      1 |      1
         |   0.00 |   0.00 |   0.00 |   0.06 |   0.06
         |   0.00 |   0.00 |   0.00 | 100.00 |
         |   0.00 |   0.00 |   0.00 |   1.00 |
---------+--------+--------+--------+--------+
Total        1618       19       22      100     1759
            91.98     1.08     1.25     5.69   100.00

*/

/* Make permenant */

data event.miso_bin_diffs_and_bayes;
  set bin_diffs_and_bayes;
run;













