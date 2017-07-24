ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Check lengths of hits against lengths of junction.
   Subset only the hits where the BLAST length is the same as the junction length, and there are no gaps
   Then count:
   novel/known * flag_junction_annotated

   Also compare against catalog! (ie, merge PB junctions to catalog junctions on coordinate,
   how many hit catalog junctions? What is the crossover with the PB type here? */


/* get junction lengths */

data junc_len;
  set evspl.splicing_events_annot_refseq;
  where flag_junction_annotated=0 and flaG_intron_retention=0;
  keep event_id event_size;
run;

data blast_hits;
  set event.unannot_junc_pacbio_hits;
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
        0 |    951 |     13 |    964
          |  37.63 |   0.51 |  38.15
          |  98.65 |   1.35 |
          |  98.65 |   0.83 |
 ---------+--------+--------+
        1 |     13 |   1550 |   1563
          |   0.51 |  61.34 |  61.85
          |   0.83 |  99.17 |
          |   1.35 |  99.17 |
 ---------+--------+--------+
 Total         964     1563     2527
             38.15    61.85   100.00

*/

proc sort data=junc_len;
  by event_id;
proc sort data=flag_top_hit;
  by event_id;
run;

data best_hits_w_junc_len;
   merge flag_top_hit (in=in1) junc_len (in=in2);
   by event_id;
   if in1 and in2;
run;


/* Check: how many hits where the alignment length is the same as the junction length?
          how many hits where the alignment length is within 10bp of junction length? 
          how many hits where the alignment length is within 20bp of junction length?*/

data flag_same_len;
  set best_hits_w_junc_len;
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
                 0         316       12.50
                 1        2211       87.50


  flag_aln_within_
            95perc    Frequency     Percent
  -------------------------------------------
                 0         209        8.27
                 1        2318       91.73


  flag_aln_within_
            90perc    Frequency     Percent
  -------------------------------------------
                 0         181        7.16
                 1        2346       92.84

Take the 2346 alignments with 90% of event length
*/

proc freq data=flag_same_len noprint;
  tables flag_aln_within_90perc*flag_best_hit*flag_first_hit / out=aln_check;
proc print data=aln_check;
run;

/*
  flag_aln_                 flag_
   within_       flag_     first_
    90perc     best_hit      hit     COUNT

       0            0          0        138
       0            0          1          1
       0            1          0          1
       0            1          1         41
       1            0          0        813
       1            0          1         12
       1            1          0         12
       1            1          1       1509

*/


/* I will carry forward the alignments where the event_size is the same as the alignment length
   and there are no gaps

   Then I will count the number of unique junctions and the total number of hits */

data subset_same_length_no_gap;
   set flag_same_len;
   if flag_aln_within_95perc=1 and flag_gapopen=0 and flag_mismatch=0 and perc_identity ge 95;
run; *2137 good alignments;

data subset_events;
  set subset_same_length_no_gap;
  keep event_id;
run;

proc sort data=subset_events nodup;
   by event_id;
run; *1399 unique events;

/* Make permenant */

data event.unannot_junc_best_blast_hits;
   set subset_same_length_no_gap;
run;

 
