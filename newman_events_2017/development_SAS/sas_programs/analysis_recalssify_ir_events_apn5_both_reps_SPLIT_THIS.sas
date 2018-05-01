ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For each cell type, flag if fragments are detected */

data set_group;
  length cell_type $3.;
  set event.mm10_refseq_intron_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type intron_id apn;
run;

data flag_intron_apn_gt0;
  set set_group;
  if apn ge 5 then flag_intron_apn_gt0=1;
  else flag_intron_apn_gt0=0;
run;

proc sort data=flag_intron_apn_gt0;
  by cell_type intron_id;
proc means data=flag_intron_apn_gt0 noprint;
  by cell_type intron_id;
  var flag_intron_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data nsc old;
  set mean_on;
  if cell_type="NSC" then output nsc;
  else output old;
run;

data flag_on_nsc;
  set nsc;
  if mean_gt0 =1 then flag_intron_nsc_on=1;
  else if mean_gt0 =0 then flag_intron_nsc_on=0;
  else flag_intron_nsc_on=.;
  keep intron_id flag_intron_nsc_on;
run;



data flag_on_old;
  set old;
  if mean_gt0 =1 then flag_intron_old_on=1;
  else if mean_gt0 =0 then flag_intron_old_on=0;
  else flag_intron_old_on=.;
  keep intron_id flag_intron_old_on;
run;

proc sort data=flag_on_nsc;
   by intron_id;
proc sort data=flag_on_old;
   by intron_id;
run;

data event.flag_intron_on_apn5;
   merge flag_on_nsc (in=in1) flag_on_old (in=in2);
   by intron_id;
   if in1 and in2;
run;


/*********************************************************/

* frag detection;


/* For each cell type, flag if fragments are detected */

data set_group;
  length cell_type $3.;
  set event.mm10_refseq_fragment_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type fragment_id apn;
run;

data flag_frag_apn_gt0;
  set set_group;
  if apn ge 5 then flag_fragment_apn_gt0=1;
  else flag_fragment_apn_gt0=0;
run;


proc sort data=flag_frag_apn_gt0;
  by cell_type fragment_id;
proc means data=flag_frag_apn_gt0 noprint;
  by cell_type fragment_id;
  var flag_fragment_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data nsc old;
  set mean_on;
  if cell_type="NSC" then output nsc;
  else output old;
run;

data flag_on_nsc;
  set nsc;
  if mean_gt0 =1 then flag_fragment_nsc_on=1;
  else if mean_gt0 =0 then flag_fragment_nsc_on=0;
  else flag_fragment_nsc_on=.;
  keep fragment_id flag_fragment_nsc_on;
run;

data flag_on_old;
  set old;
  if mean_gt0 =1 then flag_fragment_old_on=1;
  else if mean_gt0 =0 then flag_fragment_old_on=0;
  else flag_fragment_old_on=.;
  keep fragment_id flag_fragment_old_on;
run;


proc sort data=flag_on_nsc;
   by fragment_id;
proc sort data=flag_on_old;
   by fragment_id;
run;

data event.flag_fragment_on_apn5;
   merge flag_on_nsc (in=in1) flag_on_old (in=in2);
   by fragment_id;
   if in1 and in2;
run;


/*********************************************************/

* fus detection;


/* For each cell type, flag if fragments are detected */

data set_group;
  length cell_type $3.;
  set event.mm10_refseq_fusion_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type fusion_id apn;
run;

data flag_fusion_apn_gt0;
  set set_group;
  if apn ge 5 then flag_fusion_apn_gt0=1;
  else flag_fusion_apn_gt0=0;
run;


proc sort data=flag_fusion_apn_gt0;
  by cell_type fusion_id;
proc means data=flag_fusion_apn_gt0 noprint;
  by cell_type fusion_id;
  var flag_fusion_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data nsc old;
  set mean_on;
  if cell_type="NSC" then output nsc;
  else output old;
run;

data flag_on_nsc;
  set nsc;
  if mean_gt0 =1 then flag_fusion_nsc_on=1;
  else if mean_gt0 =0 then flag_fusion_nsc_on=0;
  else flag_fusion_nsc_on=.;
  keep fusion_id flag_fusion_nsc_on;
run;

data flag_on_old;
  set old;
  if mean_gt0 =1 then flag_fusion_old_on=1;
  else if mean_gt0 =0 then flag_fusion_old_on=0;
  else flag_fusion_old_on=.;
  keep fusion_id flag_fusion_old_on;
run;



proc sort data=flag_on_nsc;
   by fusion_id;
proc sort data=flag_on_old;
   by fusion_id;
run;

data event.flag_fusion_on_apn5;
   merge flag_on_nsc (in=in1) flag_on_old (in=in2);
   by fusion_id;
   if in1 and in2;
run;


/*********************************************************/



/* Merge on-flags, mean APN, introns, fragments and fusions.
   Then calculate the ratio between IR event and nearest fusion and fragment and export for plots

   I want to plot density plots of IR to fusion/fragment ratios and scatter plots of
   IR apn to fusion/fragment APN

   If IR event not detected, then we should drop the intron check */

/* Get detection flags */
data int_on;
    set event.flag_intron_on_apn5;
    where flag_intron_nsc_on ne .;
    keep intron_id flag_intron_nsc_on;
run;


data ir_on;
   set event.flag_feature_on_all_reps_ge5;
    where flag_feature_on_ge5 ne .;
   keep feature_id flag_feature_on_ge5 ;
   rename feature_id=event_id flag_feature_on_ge5=flag_splicing_on_ge5;
run;

data fus_on;
   set event.flag_fusion_on_apn5;
       where flag_fusion_nsc_on ne .;
   keep fusion_id flag_fusion_nsc_on;
run;

data frag_on;
   set event.flag_fragment_on_apn5;
       where flag_fragment_nsc_on ne .;
   keep fragment_id flag_fragment_nsc_on;
run;

/* Get mean APN */
data mean_int_apn;
   set event.mean_apn_intron_nsc;
run;

data mean_frag_apn;
   set event.mean_apn_fragment_nsc;
run;

data mean_fus_apn;
   set event.mean_apn_fusion_nsc;
run;

data mean_IR_apn;
   set event.mean_apn_IR_nsc;
run;

/* Get intron-to-IR and intron-to-fusion/fragment and merge */
data intron_to_event;
  set event.ir_events_w_intron;
   where flag_no_computed_intron=0; *only keep those with computed introns;
run;

data intron_to_fus_frag;
   set event.mm10_introns_to_fragment_fusion;
   drop chr intron_start intron_stop;
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

/* Merge all together */
proc sort data=int_on;
  by intron_id;
proc sort data=mean_int_apn;
  by intron_id;
run;
data mean_int_apn_w_dtct;
   merge int_on (in=in1) mean_int_apn (in=in2);
   by intron_id;
   if in1 and in2;
run;

proc sort data=ir_on;
  by event_id;
proc sort data=mean_ir_apn;
  by event_id;
run;
data mean_ir_apn_w_dtct;
   merge ir_on (in=in1) mean_ir_apn (in=in2);
   by event_id;
   if in1 and in2;
run;

proc sort data=fus_on;
  by fusion_id;
proc sort data=mean_fus_apn;
  by fusion_id;
run;
data mean_fus_apn_w_dtct;
   merge fus_on (in=in1) mean_fus_apn (in=in2);
   by fusion_id;
   if in1 and in2;
run;


proc sort data=frag_on;
  by fragment_id;
proc sort data=mean_frag_apn;
  by fragment_id;
run;
data mean_frag_apn_w_dtct;
   merge frag_on (in=in1) mean_frag_apn (in=in2);
   by fragment_id;
   if in1 and in2;
run;

data frag_5prime_mean_dtct;
   set mean_frag_apn_w_dtct;
   rename fragment_id=fragment_id_5prime flag_fragment_nsc_on=flag_fragment_on_5prime
          mean_apn_fragment=mean_apn_fragment_5prime;
run;

data frag_3prime_mean_dtct;
   set mean_frag_apn_w_dtct;
   rename fragment_id=fragment_id_3prime flag_fragment_nsc_on=flag_fragment_on_3prime
          mean_apn_fragment=mean_apn_fragment_3prime;
run;

data fus_5prime_mean_dtct;
   set mean_fus_apn_w_dtct;
   rename fusion_id=fusion_id_5prime flag_fusion_nsc_on=flag_fusion_on_5prime
          mean_apn_fusion=mean_apn_fusion_5prime;
run;

data fus_3prime_mean_dtct;
   set mean_fus_apn_w_dtct;
   rename fusion_id=fusion_id_3prime flag_fusion_nsc_on=flag_fusion_on_3prime
          mean_apn_fusion=mean_apn_fusion_3prime;
run;


proc sort data=frag_5prime_mean_dtct;
   by fragment_id_5prime;
proc sort data=frag_3prime_mean_dtct;
   by fragment_id_3prime;
proc sort data=fus_5prime_mean_dtct;
   by fusion_id_5prime;
proc sort data=fus_3prime_mean_dtct;
   by fusion_id_3prime;
proc sort data=mean_ir_apn_w_dtct;
   by event_id;
proc sort data=mean_int_apn_w_dtct;
   by intron_id;
proc sort data=intron_to_ir_fus_frag;
   by intron_id;
run;

data intron2ir_p1;
  merge intron_to_ir_fus_frag (in=in1) mean_int_apn_w_dtct (in=in2);
  by intron_id;
  if in1 and in2;
run;

proc sort data=intron2ir_p1;
   by event_id;
run;

data intron2ir_p2;
   merge intron2ir_p1 (in=in1) mean_ir_apn_w_dtct (in=in2);
   by event_id;
   if in1 and in2;
run;

proc sort data=intron2ir_p2;
  by fusion_id_5prime;
run;

data intron2ir_p3;
  merge intron2ir_p2 (in=in1) fus_5prime_mean_dtct (in=in2);
  by fusion_id_5prime;
  if in1 and in2;
run;

proc sort data=intron2ir_p3;
  by fusion_id_3prime;
run;


data intron2ir_p4;
  merge intron2ir_p3 (in=in1) fus_3prime_mean_dtct (in=in2);
  by fusion_id_3prime;
  if in1 and in2;
run;

proc sort data=intron2ir_p4;
  by fragment_id_5prime;
run;

data intron2ir_p5;
  merge intron2ir_p4 (in=in1) frag_5prime_mean_dtct (in=in2);
  by fragment_id_5prime;
  if in1 and in2;
run;

proc sort data=intron2ir_p5;
  by fragment_id_3prime;
run;

data intron2ir_p6;
  merge intron2ir_p5 (in=in1) frag_3prime_mean_dtct (in=in2);
  by fragment_id_3prime;
  if in1 and in2;
run;

/* Calc IR ratios */
data ir_ratios;
  set intron2ir_p6;
  ir_to_5fus_ratio=mean_apn_ir/mean_apn_fusion_5prime;
  ir_to_5frag_ratio=mean_apn_ir/mean_apn_fragment_5prime;
  ir_to_3fus_ratio=mean_apn_ir/mean_apn_fusion_3prime;
  ir_to_3frag_ratio=mean_apn_ir/mean_apn_fragment_3prime;
  mean_apn_fusion_mean=(mean_apn_fusion_5prime+mean_apn_fusion_3prime)/2;
  mean_apn_fragment_mean=(mean_apn_fragment_5prime+mean_apn_fragment_3prime)/2;
  ir_to_mean_fus_ratio=mean_apn_ir/mean_apn_fusion_mean;
  ir_to_mean_frag_ratio=mean_apn_ir/mean_apn_fragment_mean;
run;

/* Make permenant */

data event.intron2ir_frag_fus_all_info_apn5;
   set ir_ratios;
run;



/* Flag if in expressed gene and export */

data exp_gene;
  set event.flag_gene_expressed;
  keep gene_id flag_gene_expressed;
run;

proc sort data=exp_gene;
  by gene_id;
proc sort data=ir_ratios;
  by gene_id;
run;

data ir_ratios_w_gene;
  merge ir_ratios (in=in1) exp_gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc export data=ir_ratios_w_gene outfile="!MCLAB/event_analysis/analysis_output/event_analysis_intron_ir_analysis2.csv"
   dbms=csv replace;
run;



/*********************************************************/

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
   set event.intron2ir_frag_fus_all_info_apn5;
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
        mean_apn_intron flag_splicing_on_ge5 mean_apn_ir flag_fusion_on mean_apn_fusion flag_fragment_on
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
   if flag_fusion_on=1 and flag_splicing_on_ge5=1 then flag_ir_event_analyzable=1;
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

/*********************************************************/


/* Reclassifying IR events.

I am going to try the following criteria:
(1) If IR >~ donor exon (90%), then novel donor
(2) Else if IR > ~10 donor then IR
(3) Else unprocessed

Do for IR and for intron and check the concordance */

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
   set event.intron2ir_frag_fus_all_info_apn5;
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
        mean_apn_intron flag_splicing_on_ge5 mean_apn_ir flag_fusion_on mean_apn_fusion flag_fragment_on
        mean_apn_fragment ir_to_fus_ratio ir_to_frag_ratio;
run;


/*

(1) If IR APN > 5 and IR  ~ exon, then flag_possible_novel_donor=1
(2) If IR APN > 5, IR ~ intron, then flag_possible_IR
(3) If IR APN > 5, IR > intron, then flag_unprocessed_rna=1
(4) If IR APN < 5 then flag_low_expression=1

Figure 1. Possible classifications for putative IR events. Putative IR events are considered representative of a possible novel donor site if the abundance of the event is equivalent to the abundance of its adjacent 5’ exonic region and there is little expression of the adjacent intron. Events are considered as possible IR events if the abundance of the putative IR event and the adjacent intron are equivalent. The remaining set of putative IR events consist of instances where there is low abundance of the event and the intron relative to the 5’ exon, and are considered possible unprocessed transcript.

*/



data classify_introns_ir2;
  set intron2ir_min_info;
   if flag_fusion_on=1 and flag_splicing_on_ge5=1 then do; *only do for "analyzable" events;
   flag_low_expressed=0;
       /* APN > 1*/
     if mean_apn_ir ge 5 then do;
       flag_possible_unprocessed=0;
       if mean_apn_ir ge (0.9 * mean_apn_fusion) then flag_possible_novel_donor=1;
       else flag_possible_novel_donor=0;
       if mean_apn_ir ge (0.1 * mean_apn_fusion) then do;
              if mean_apn_intron ge (0.5 * mean_apn_IR) then do;
                    flag_possible_ir=1;
                    flag_possible_unprocessed=0; end;
              else do;
                    flag_possible_ir=0;
                    flag_possible_unprocessed=1; end;
        end;
        else do;
          flag_possible_ir=0;
          flag_possible_novel_donor=0;
          flag_possible_unprocessed=1; end;
         end;
     else  do;
          flag_possible_ir=0;
          flag_possible_novel_donor=0;flag_possible_unprocessed=1;
     end;  end;
   else flag_low_expressed=1;
run;

proc freq data=classify_introns_ir2 noprint;
   where flag_fusion_on=1 and flag_splicing_on_ge5=1;
   tables flag_low_expressed*flag_possible_novel_donor*flag_possible_ir*flag_possible_unprocessed / out=ir_class;
run;

proc print data=ir_class;
run;

data check;
   set classify_introns_ir2;
   where flag_possible_novel_donor=0 and flag_possible_ir=0 and flag_possible_unprocessed=0;
run;



/*
                                    flag_
   flag_low_    flag_possible_    possible_    flag_possible_
   expressed      novel_donor         ir         unprocessed     COUNT    PERCENT

       0               0              0               1          24794    97.7990
       0               0              1               0            200     0.7889
       0               1              0               1            314     1.2386
       0               1              1               0             44     0.1736



*/



proc freq data=classify_introns_ir2;
   tables flag_splicing_on_ge5*flag_fusion_on;
run;

*1669 events where the donor fusion is off but the IR is on;
*25352 events where the donor fusion and IR are on;


/* Export APNs and flags so I can make plots of IR vs donor, IR vs intron, intron vs donor */

data data_for_plots;
  set classify_introns_ir2;
  keep event_id flag_no_computed_intron mean_apn_intron flag_splicing_on_ge5 mean_apn_ir
       flag_fusion_on mean_apn_fusion flag_low_expressed flag_possible_novel_donor flag_possible_ir flag_possible_unprocessed ;
run;

proc export data=data_for_plots
    outfile="!MCLAB/event_analysis/analysis_output/event_analysis_reclassified_ir_jrbn3_apn5_both_reps.csv"
   dbms=csv replace;
run;



/* Make permenant */

data event.ir_reclassification_v3_apn5;
   set classify_introns_ir2;
run;


/*********************************************************/


/* Count reclassified IR events by type and by whether they have a PB hit or not */

* IR with hit;
data ir_w_hit;
  set event.blast_junc_to_pb_w_annot;
  where flag_intron_retention=1 and flag_feature_on_ge5=1
        and flag_event_has_hit=1;
  keep event_id;
run;

* IR reclassification;
data ir_class;
  set event.ir_reclassification_v3_apn5;
  where flag_splicing_on_ge5=1 and flag_low_expressed=0;
  keep event_id flag_low_expressed flag_possible_unprocessed flag_possible_novel_donor flag_possible_ir ;
run;


proc sort data=ir_w_hit nodup;
   by event_id;
proc sort data=ir_class;
   by event_id;
run;

data ir_class_pb_hit;
  merge ir_class (in=in1) ir_w_hit (in=in2);
  by event_id;
  if in2 then flag_pb_hit=1; else flag_pb_hit=0;
  if in1 then output;
run;


proc freq data=ir_class_pb_hit noprint;
   tables flag_low_expressed*flag_possible_unprocessed*
          flag_possible_novel_donor*flag_possible_ir*flag_pb_hit / out=ir_count;
run;

proc print data=ir_count;
run;


/*
                                                    flag_
 flag_low_    flag_possible_    flag_possible_    possible_    flag_pb_
 expressed      unprocessed       novel_donor         ir          hit      COUNT

     0               0                 0              1            0         71
     0               0                 0              1            1         52
     0               0                 1              1            0          5
     0               0                 1              1            1          6
     0               1                 0              0            0        240
     0               1                 0              0            1         94
     0               1                 1              0            0         29
     0               1                 1              0            1          8

TABLE

CLASS           NO HIT  HIT
Poss. IR        71      52
Ambig IR        5       6
Unprocessed     240     94
Novel donor     29      8

*/


