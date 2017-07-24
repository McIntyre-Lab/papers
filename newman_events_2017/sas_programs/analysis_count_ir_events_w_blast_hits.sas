ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* If an IR event has at least one hit, then flag the IR event has having Pacbio support
   and look at the cross-tabulation with IR classification

   Then count the number of PB transcripts by type for each IR event */

data ir_hits;
 set event.unannot_ir_best_blast_hits;
 if perc_identity < 0.95 then delete;
run;


data ir_hits_uniq;
  set ir_hits;
  keep event_id flag_possible_unprocessed flag_possible_novel_donor flag_possible_ir ;
run;

proc sort data=ir_hits_uniq nodup;
   by event_id;
proc freq data=ir_hits_uniq noprint;
   tables flag_possible_unprocessed*flag_possible_novel_donor*flag_possible_ir / out=ir_counts;
proc print data=ir_counts;
run;

/*
                                        flag_
  flag_possible_    flag_possible_    possible_
    unprocessed       novel_donor         ir       COUNT    NOTE

         0                 0              1           90    Possible IR
         0                 1              1           25    Ambiguous IR
         1                 0              0         1052    Unprocessed RNA
         1                 1              0           42    Possible novel donor

*/

/* Merge with full set of classified IR events and look at cross-over */

data ir_blasted;
   set event.ir_event_pacbio_hits;
   keep event_id;
run;

data ir_class;
  set event.ir_reclassification_v2;
  keep event_id flag_low_expressed flag_possible_unprocessed flag_possible_novel_donor flag_possible_ir;
run;

data ir_hits_uniq2;
  set ir_hits_uniq;
  keep event_id;
run;

proc sort data=ir_class nodup;
  by event_id;
proc sort data=ir_blasted nodup;
  by event_id;
proc sort data=ir_hits_uniq2 nodup;
  by event_id;
run;

data ir_class_flag_hits;
   merge ir_class (in=in1) ir_blasted (in=in2) ir_hits_uniq2 (in=in3);
   by event_id;
   if in2 then flag_ir_blast_hit=1; else flag_ir_blast_hit=0;
   if in3 then flag_ir_good_hit=1; else flag_ir_good_hit=0;
   if in1 then output;
run;

proc freq data=ir_class_flag_hits noprint;
   tables flag_low_expressed*flag_possible_unprocessed*flag_possible_novel_donor* flag_possible_ir*flag_ir_blast_hit*flag_ir_good_hit / out=ir_count_by_class;
proc print data=ir_count_by_class;
run;


/*
 flag_low_  flag_possible_  flag_possible_  possible_   flag_ir_  flag_ir_
 expressed    unprocessed     novel_donor       ir     blast_hit  good_hit   COUNT  PERCENT

     0             0               0            1          0          0        180   0.5175
     0             0               0            1          1          0          3   0.0086
     0             0               0            1          1          1         90   0.2588
     0             0               1            1          0          0         30   0.0863
     0             0               1            1          1          1         25   0.0719
     0             1               0            0          0          0      32854  94.4623
     0             1               0            0          1          0        161   0.4629
     0             1               0            0          1          1       1052   3.0247
     0             1               1            0          0          0        254   0.7303
     0             1               1            0          1          0         89   0.2559
     0             1               1            0          1          1         42   0.1208
     1             .               .            .          0          0     193187    .



SUMMARY:
Type		NO BLAST	NO GOOD HIT	GOOD HIT
Low expression	193187		n/a		n/a		(note: these were not blasted)
Unprocessed	32854		161		1052
Novel Donor	254		89		42
Poss. IR	180		3		90
IR/Novel donor	30		0		25

*/

/* Now I want to count by PB transcript type
   Basically IR class * PB types

   So, for each IR, take the list of "good" hits and count the number of PB transcript types
   Then crosstabs on IR type and PB type list
   */

proc sort data=ir_hits;
  by event_id;
proc sort data=ir_class;
  by event_id;
run;

data ir_hits_w_class;
  merge ir_hits (in=in1) ir_class (in=in2);
  by event_id;
  if in1 and in2;
run;

data ir_hits2class;
  set ir_hits_w_class;
  keep event_id pb_status;
run;

proc sort data=ir_hits2class nodup;
   by event_id pb_status;
proc freq data=ir_hits2class noprint;
   tables event_id / out=pb_count;
proc sort data=pb_count;
  by descending count;
run; *max 2;


data cat_pb; 
  array pb[2] $ 30;
  retain pb1-pb2;
  set ir_hits2class;
  by event_id;
  if first.event_id then do;
     call missing(of pb1-pb2);
     records = 0;
  end;
  records + 1;
  pb[records]=pb_status;
  if last.event_id then output;
run;

  *clean up the output file;
data cat_pb2;
  set cat_pb;
  length pb_hit_type $ 100;
         pb_hit_type= catx("|", OF pb1-pb2);
  keep event_id pb_hit_type;
  run;

/* Merge in IR class flags */

data ir_class_w_cat_pb;
  merge cat_pb2 (in=in1) ir_class (in=in2);
  by event_id;
  if in1;
run;

data ir_class_w_cat_pb2;
  set ir_class_w_cat_pb;
  length ir_class $50.;
  if flag_low_expressed=1 then ir_class="Low expression";
  else if flag_possible_novel_donor=1 and flag_possible_ir=1 then ir_class="Ambiguous IR";
  else if flag_possible_novel_donor=0 and flag_possible_ir=1 then ir_class="Possible IR";
  else if flag_possible_novel_donor=1 and flag_possible_ir=0 then ir_class="Novel donor";
  else ir_class="Unprocessed";
run;



proc freq data=ir_class_w_cat_pb2;
   tables ir_class*pb_hit_type;
run;


/*
IE, where are the best hits ?

ir_class      pb_hit_type

Frequency    |
Percent      |
Row Pct      |
Col Pct      |Known   |Known|No|Novel   |  Total
             |        |vel     |        |
-------------+--------+--------+--------+
Ambiguous IR |      5 |      6 |     14 |     25
             |   0.41 |   0.50 |   1.16 |   2.07
             |  20.00 |  24.00 |  56.00 |
             |   2.40 |   8.00 |   1.51 |
-------------+--------+--------+--------+
Novel donor  |      9 |      8 |     25 |     42
             |   0.74 |   0.66 |   2.07 |   3.47
             |  21.43 |  19.05 |  59.52 |
             |   4.33 |  10.67 |   2.70 |
-------------+--------+--------+--------+
Possible IR  |     17 |     17 |     56 |     90
             |   1.41 |   1.41 |   4.63 |   7.44
             |  18.89 |  18.89 |  62.22 |
             |   8.17 |  22.67 |   6.05 |
-------------+--------+--------+--------+
Unprocessed  |    177 |     44 |    831 |   1052
             |  14.64 |   3.64 |  68.73 |  87.01
             |  16.83 |   4.18 |  78.99 |
             |  85.10 |  58.67 |  89.74 |
-------------+--------+--------+--------+
Total             208       75      926     1209
                17.20     6.20    76.59   100.00

Most hits are to novel transcripts, so this is okay

*/

