ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Setting a bunch of flags so I can explore the IR/intron data

   (1) If the donor exonic region of an IR event is not detected, then flag IR as not analyzable
       Else, IR is analyzable
   (2) Drop IR events that are not detected
   (3) Bin IR events into APN <1, APN <5, APN<10 and APN>10, and also log(IR-apn/exon-apn) <-5 and >-5
       and check cross-tabulation. This should give me an idea as to where to set a coverage threshold
       for determining if an IR event represents unprocessed transcript or not
   (4) Flag if intron APN > 1, >5, >10. I want to use this as evidence of intron retention
       vs evidence of possible novel donor.

*/

/* Get intron-to-IR and intron-to-fusion/fragment and merge */
data intron_to_event;
  set event.ir_events_w_intron;
   where flag_no_computed_intron=0; *only keep those with computed introns;
run;

data intron_to_fus_frag;
   set event.mm10_introns_to_fragment_fusion;
run;

data ir_coord;
   set evspl.splicing_events_annot_refseq;
   where flag_intron_retention=1;
   keep event_id chr feature1_start feature1_stop
        feature2_start feature2_stop;
run;

proc sort data=intron_to_event;
   by intron_id;
proc sort data=intron_to_fus_frag;
   by intron_id;
run;

data intron_to_ir_fus_frag;
  merge intron_to_event (in=in1) intron_to_fus_frag (in=in2);
  by intron_id;
  if in1 and in2;
run;

proc sort data=intron_to_ir_fus_frag;
  by event_id;
proc sort data=ir_coord;
  by event_id;
run;

data intron_to_ir_W_coord;
  merge intron_to_ir_fus_frag (in=in1) ir_coord (in=in2);
  by event_id;
  if in1 and in2;
  if intron_start = feature2_start-1 then flag_ir_from_3prime=1; else flag_ir_from_3prime=0;
  if intron_stop = feature1_stop-1 then flag_ir_from_5prime=1; else flag_ir_from_5prime=0;
  keep event_id flag_ir_from_3prime flag_ir_from_5prime;
run;

proc freq data=intron_to_ir_w_coord;
   tables flag_ir_from_5prime*flag_ir_from_3prime;
run;

/* First, I want keep only the corresponding donor fusion/fragment for each IR event */

data intron2ir_info;
   set event.intron2ir_frag_fus_all_info;
run;

proc sort data=intron2ir_info;
  by event_id;
proc sort data=intron_to_ir_w_coord;
  by event_id;
run;

data intron2ir_info2;
  merge intron2ir_info (in=in1) intron_to_ir_w_coord (in=in2);
  by event_id;
  if in1 and in2;
run;

data intron2ir_min_info;
   set intron2ir_info2;
   length fusion_id $10.;
   length fragment_id $20.;
   if flag_ir_from_5prime=1 then do;
     fusion_id=fusion_id_5prime;
     fragment_id=fragment_id_5prime;
     flag_fusion_on=flag_fusion_on_5prime;
     mean_apn_fusion=mean_apn_fusion_5prime;
     flag_fragment_on=flag_fragment_on_5prime;
     mean_apn_fragment=mean_apn_fragment_5prime;
     ir_to_fus_ratio=ir_to_5fus_ratio;
     ir_to_frag_ratio=ir_to_5frag_ratio; end;
   else do;
     fusion_id=fusion_id_3prime;
     fragment_id=fragment_id_3prime;
     flag_fusion_on=flag_fusion_on_3prime;
     mean_apn_fusion=mean_apn_fusion_3prime;
     flag_fragment_on=flag_fragment_on_3prime;
     mean_apn_fragment=mean_apn_fragment_3prime;
     ir_to_fus_ratio=ir_to_3fus_ratio;
     ir_to_frag_ratio=ir_to_3frag_ratio; end;
   keep event_id intron_id gene_id flag_no_computed_intron fusion_id fragment_id flag_intron_nsc_on
        mean_apn_intron flag_splicing_on mean_apn_ir flag_fusion_on mean_apn_fusion flag_fragment_on
        mean_apn_fragment ir_to_fus_ratio ir_to_frag_ratio;
run;

/* Set flags:
I am going to try the following criteria:
(1) If IR >~ donor exon (90%), then novel donor
(2) Else if IR > ~10 donor then IR
(3) Else unprocessed

Do for IR and for intron and check the concordance
*/

/* calc median apn ratio */

proc means data=intron2ir_min_info;
   where ir_to_fus_ratio > 0;
   var ir_to_fus_ratio;
   output out=ir_median_ratio median=;
run;

proc print data=ir_median_ratio;
run; *median ratio = 0.028751;


data set_ir_flags;
   set intron2ir_min_info;
   /* IF fusion and IR event are detected, then event is analyzable */
   if flag_fusion_on=1 and flag_splicing_on=1 then flag_ir_event_analyzable=1;
   else flag_ir_event_analyzable=0;

   /* Set IR APN bins */
   if mean_apn_ir < 1 then ir_apn_bin=1;
   else if mean_apn_ir < 5 then ir_apn_bin=2;
   else if mean_apn_ir < 10 then ir_apn_bin=3;
   else if mean_apn_ir ge 10 then ir_apn_bin=4;

   /* Set IR-to-exon ration flag */
   if ir_to_fus_ratio ge 0.029521 then flag_ir_ratio_ge_median=1;
   else flag_ir_ratio_ge_median=0;

   /* Set intron APN bins */
   if mean_apn_intron < 1 then intron_apn_bin=1;
   else if mean_apn_intron < 5 then intron_apn_bin=2;
   else if mean_apn_intron < 10 then intron_apn_bin=3;
   else if mean_apn_intron ge 10 then intron_apn_bin=4;
run;

/* CHeck freqs */

proc freq data=set_ir_flags;
  tables flag_ir_event_analyzable;
run;

* 25352 of 154934 events analyzable;

proc freq data=set_ir_flags;
  where flag_ir_event_analyzable=1;
  tables flag_ir_ratio_ge_median*ir_apn_bin
         flag_ir_ratio_ge_median*intron_apn_bin
         ir_apn_bin*intron_apn_bin;
run;

/*          Table of flag_ir_ratio_ge_median by ir_apn_bin

 flag_ir_ratio_ge_median     ir_apn_bin

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|       2|       3|       4|  Total
 ---------+--------+--------+--------+--------+
        0 |  10064 |   2477 |    190 |     83 |  12814
          |  39.70 |   9.77 |   0.75 |   0.33 |  50.54
          |  78.54 |  19.33 |   1.48 |   0.65 |
          |  53.53 |  46.02 |  31.05 |  14.95 |
 ---------+--------+--------+--------+--------+
        1 |   8738 |   2906 |    422 |    472 |  12538
          |  34.47 |  11.46 |   1.66 |   1.86 |  49.46
          |  69.69 |  23.18 |   3.37 |   3.76 |
          |  46.47 |  53.98 |  68.95 |  85.05 |
 ---------+--------+--------+--------+--------+
 Total       18802     5383      612      555    25352
             74.16    21.23     2.41     2.19   100.00

flag_ir_ratio_ge_median     intron_apn_bin

Frequency|
Percent  |
Row Pct  |
Col Pct  |       1|       2|       3|       4|  Total
---------+--------+--------+--------+--------+
       0 |   7382 |   4315 |    682 |    435 |  12814
         |  29.12 |  17.02 |   2.69 |   1.72 |  50.54
         |  57.61 |  33.67 |   5.32 |   3.39 |
         |  46.16 |  59.64 |  55.36 |  48.71 |
---------+--------+--------+--------+--------+
       1 |   8610 |   2920 |    550 |    458 |  12538
         |  33.96 |  11.52 |   2.17 |   1.81 |  49.46
         |  68.67 |  23.29 |   4.39 |   3.65 |
         |  53.84 |  40.36 |  44.64 |  51.29 |
---------+--------+--------+--------+--------+
Total       15992     7235     1232      893    25352
            63.08    28.54     4.86     3.52   100.00

        Table of ir_apn_bin by intron_apn_bin

ir_apn_bin     intron_apn_bin

Frequency|
Percent  |
Row Pct  |
Col Pct  |       1|       2|       3|       4|  Total
---------+--------+--------+--------+--------+
       1 |  13942 |   4442 |    276 |    142 |  18802
         |  54.99 |  17.52 |   1.09 |   0.56 |  74.16
         |  74.15 |  23.63 |   1.47 |   0.76 |
         |  87.18 |  61.40 |  22.40 |  15.90 |
---------+--------+--------+--------+--------+
       2 |   1697 |   2580 |    790 |    316 |   5383
         |   6.69 |  10.18 |   3.12 |   1.25 |  21.23
         |  31.53 |  47.93 |  14.68 |   5.87 |
         |  10.61 |  35.66 |  64.12 |  35.39 |
---------+--------+--------+--------+--------+
       3 |    131 |    125 |    129 |    227 |    612
         |   0.52 |   0.49 |   0.51 |   0.90 |   2.41
         |  21.41 |  20.42 |  21.08 |  37.09 |
         |   0.82 |   1.73 |  10.47 |  25.42 |
---------+--------+--------+--------+--------+
       4 |    222 |     88 |     37 |    208 |    555
         |   0.88 |   0.35 |   0.15 |   0.82 |   2.19
         |  40.00 |  15.86 |   6.67 |  37.48 |
         |   1.39 |   1.22 |   3.00 |  23.29 |
---------+--------+--------+--------+--------+
Total       15992     7235     1232      893    25352
            63.08    28.54     4.86     3.52   100.00

Proposed IR re-classification:

(a) If IR_APN > 1 and IR_APN_ratio > median then event is either IR, or is novel exon
(b) If IR_APN < 1 and IR_APN_ratio > median then unresolvable (ie, expression is too low to tell)
(c) If IR_APN > 1 and IR_APN_ratio < median then unprocessed transcript
(d) If IR_APN < 1 and IR_APN_ratio < median then unresolvable (ie, expression is too low to tell)

Then for (a), if intron APN > 1 then likely IR
Else likely novel donor

*/
