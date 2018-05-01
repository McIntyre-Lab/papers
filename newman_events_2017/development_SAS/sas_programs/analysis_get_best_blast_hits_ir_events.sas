ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Need to match these up to reclassified IR types and check lengths of hits against lengths of event.
   Subset only the hits where the BLAST length is the same as the junction length, and there are no gaps
   Then count:
   novel/known * ir type */


/* get junction lengths */

data ir_len;
  set evspl.splicing_events_annot_refseq;
  where flaG_intron_retention=1;
  keep event_id event_size;
run;

data ir_class;
   set event.ir_reclassification_v2;
   where flag_low_expressed = 0;
   keep event_id flag_possible_unprocessed flag_possible_novel_donor flag_possible_ir;
run;

proc sort data=ir_len;
  by event_id;
proc sort data=ir_class;
  by event_id;
run;

data ir_len2;
  merge ir_len (in=in1) ir_class (in=in2);
  by event_id;
  if in1 and in2;
run;



data blast_hits;
  set event.ir_event_pacbio_hits;
run;

proc sort data=blast_hits;
   by event_id;
run;

data blast_hits2;
  set blast_hits;
  by event_id;
  if first.event_id then flag_first_hit=1;
  else flag_first_hit=0;
run;

/* Sort by length, e-value and identity, then flag_best_hit */

proc sort data=blast_hits2;
   by event_id descending length evalue descending perc_identity;
run;

data flag_top_hit;
  set blast_hits2;
  by event_id;
  if first.event_id then flag_best_hit=1;
  else flag_best_hit=0;
run;

proc freq data=flag_top_hit;
   tables flag_best_hit*flag_first_hit;
run;

/*

flag_best_hit     flag_first_hit

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   3311 |     58 |   3369
         |  68.54 |   1.20 |  69.74
         |  98.28 |   1.72 |
         |  98.28 |   3.97 |
---------+--------+--------+
       1 |     58 |   1404 |   1462
         |   1.20 |  29.06 |  30.26
         |   3.97 |  96.03 |
         |   1.72 |  96.03 |
---------+--------+--------+
Total        3369     1462     4831
            69.74    30.26   100.00

*/

proc sort data=ir_len2;
  by event_id;
proc sort data=flag_top_hit;
  by event_id;
run;

data best_hits_w_ir_len;
   merge flag_top_hit (in=in1) ir_len2 (in=in2);
   by event_id;
   if in1 and in2;
run;


/* Check: how many hits where the alignment length is the same as the junction length?
          how many hits where the alignment length is within 95% of junction length? 
          how many hits where the alignment length is within 90% of junction length?*/

data flag_same_len;
  set best_hits_w_ir_len;
  if length=event_size then flag_aln_eq_event_size=1;
  else flag_aln_eq_event_size=0;
  if length ge event_size*0.95 then flag_aln_within_95perc=1;
  else flag_aln_within_95perc=0;
  if length ge event_size*0.9 then flag_aln_within_90perc=1;
  else flag_aln_within_90perc=0;
run;

proc freq data=flag_same_len;
  tables flag_aln_eq_event_size flag_aln_within_95perc flag_aln_within_90perc;
run;


/*
     flag_aln_eq_
       event_size    Frequency     Percent
 -------------------------------------------
                0        4682       96.92
                1         149        3.08


 flag_aln_within_
           95perc    Frequency     Percent
 -------------------------------------------
                0        2068       42.81
                1        2763       57.19


 flag_aln_within_
           90perc    Frequency     Percent
 -------------------------------------------
                0        1695       35.09
                1        3136       64.91


Take the 3136 alignments 90% of the event size
*/

proc freq data=flag_same_len noprint;
  tables flag_aln_within_90perc*flag_best_hit*flag_first_hit / out=aln_check;
proc print data=aln_check;
run;

/*
 flag_aln_                 flag_
  within_       flag_     first_
   90perc     best_hit      hit     COUNT    PERCENT

     0            0          0       1613    33.3885
     0            0          1         14     0.2898
     0            1          0         12     0.2484
     0            1          1         56     1.1592
     1            0          0       1698    35.1480
     1            0          1         44     0.9108
     1            1          0         46     0.9522
     1            1          1       1348    27.9031
*/


/* I will carry forward the alignments where the alignment length is 95% of the event szie
   and there are no gaps

   Then I will count the number of unique junctions and the total number of hits */

data subset_same_length_no_gap;
   set flag_same_len;
   if flag_aln_within_90perc=1 and flag_gapopen=0 and flag_mismatch=0 and perc_identity ge 95;
run; *1580 good alignments;

data subset_events;
  set subset_same_length_no_gap;
  keep event_id;
run;

proc sort data=subset_events nodup;
   by event_id;
run; *1209 unique events;

/* Make permenant */

data event.unannot_ir_best_blast_hits;
   set subset_same_length_no_gap;
run;

 

