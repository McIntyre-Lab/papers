/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Simulated transcripts vs Refseq BLAST: Identify what transcripts were simulated if possible

   BLAST processing outline:
   if query has multiple hits to the same transcript, then remove -> fragmented hit
   if query length = ref length, no mismatches or gaps -> "True hit"
   query = ref length, with no more than 5 mismatches -> "True hit with mismatches"
   diff(query, ref)  < 5bp, no mismatches -> "Good hit, 5bp diff"
   diff(query, ref)  < 5bp, no more than 5 mismatches -> "Good hit, 5bp diff, mismatches"
   diff(query, ref)  < 10bp, no mismatches -> "Good hit, 10bp diff"
   diff(query, ref)  < 10bp, no more than 5 mismatches -> "Good hit, 10bp diff, mismatches"
   diff(query, ref) ~ 90%, no mismatches -> "90% hit"
   diff(query, ref) ~ 90%, no more than 5 mismatches -> "90% hit, mismatches"
   Take the best bitscore as the final

   I am not allowing for any gaps, but may allow for up to 5 mismatches

   When EACH of these is met, remove all other hits for that query transcript


   Macrotize!

*/

/* Identify fragmented hits */

data blast_hits;
  set event.sim_xscript2refseq_blast;
run;

proc sort data=blast_hits;
   by simulation test query_id transcript_id;
proc freq data=blast_hits noprint;
   by simulation test query_id;
   tables transcript_id / out=ref_hits_per_query;
run;

data fragmented;
  set ref_hits_per_query;
  where count > 1;
  keep simulation test query_id transcript_id;
run;

proc sort data=blast_hits;
   by simulation test query_id transcript_id;
proc sort data=fragmented;
   by simulation test query_id transcript_id;
run;

data blast_hits_no_frag;
  merge blast_hits (in=in1) fragmented (in=in2);
  by simulation test query_id transcript_id;
  if in1 and in2 then delete;
  else if in1 then output;
run; *1369725 blast hits in total, reduced to 493656;

/* Set up flags */

data match_flags;
   set blast_hits_no_frag;
   /* Flag mismatches and gaps first */
   if n_mismatch > 0 then flag_has_mismatch=1; else flag_has_mismatch=0;
   if n_mismatch > 5 then flag_mismatch_gt5=1; else flag_mismatch_gt5=0;
   if n_gapopen > 0 then flag_gapopen=1; else flag_gapopen=0;

   /* True match -- both query and reference have the same length */
   if query_seq_length = transcript_seq_length and flag_has_mismatch=0
   then flag_len_query_eq_ref=1;
   else flag_len_query_eq_ref=0;

   /* True match - with mismatches */
   if query_seq_length = transcript_seq_length and flag_has_mismatch=1
   and flag_mismatch_gt5=0 then flag_len_query_eq_ref_mm5=1;
   else flag_len_query_eq_ref_mm5=0;

   /* Complete match -- BLAST encompasses entire reference sequence*/
   if hit_length = transcript_seq_length and flag_has_mismatch=0
   then flag_len_match_eq_ref=1;
   else flag_len_match_eq_ref=0;

   /* Complete match with mismatches */
   if hit_length = transcript_seq_length and flag_has_mismatch=1
   and flag_mismatch_gt5=0 then flag_len_match_eq_ref_mm5=1;
   else flag_len_match_eq_ref_mm5=0;

   /* match and ref differ by 5 bp */
   if abs(hit_length-transcript_seq_length) le 5  and flag_has_mismatch=0 
   then flag_len_match_5bp_diff_ref=1;
   else flag_len_match_5bp_diff_ref=0;

   /* match and ref differ by 5 bp -- with mismatches*/
   if abs(hit_length-transcript_seq_length) le 5  and flag_has_mismatch=1 
   and flag_mismatch_gt5=0 then flag_len_match_5bp_diff_ref_mm5=1;
   else flag_len_match_5bp_diff_ref_mm5=0;

   /* match and ref differ by 10 bp */
   if abs(hit_length-transcript_seq_length) le 10  and flag_has_mismatch=0 
   then flag_len_match_10bp_diff_ref=1;
   else flag_len_match_10bp_diff_ref=0;

   /* match and ref differ by 10 bp -- with mismatches */
   if abs(hit_length-transcript_seq_length) le 10  and flag_has_mismatch=1 
   and flag_mismatch_gt5=0 then flag_len_match_10bp_diff_ref_mm5=1;
   else flag_len_match_10bp_diff_ref_mm5=0;

   /* match and ref within 90% of total length of ref */
   if hit_length/transcript_seq_length ge 0.9 and flag_has_mismatch = 0
   then flag_len_match_90prc_ref=1;
   else flag_len_match_90prc_ref=0;

   /* match and ref within 90% of total length of ref -- with mismatches */
   if hit_length/transcript_seq_length ge 0.9 and flag_has_mismatch = 1
   and flag_mismatch_gt5=0 then flag_len_match_90prc_ref_mm5=1;
   else flag_len_match_90prc_ref_mm5=0;
run;

/* Macro to pull out hits of a given type, and remove all other hits of that query */

%macro extractHits(datain,flag,hitout,otherout);

data hits others;
   set &datain.;
   if flag_gapopen=0 and &flag.=1 then output hits;
   else output others;
run;

data query2drop;
  set hits;
  keep simulation test query_id;
run;

proc sort data=others;
  by simulation test query_id;
proc sort data=query2drop nodup;
  by simulation test query_id;
run;

data &otherout.;
  merge others (in=in1) query2drop (in=in2);
  by simulation test query_id;
  if in1 and in2 then delete;
  else if in1 then output;
run;

data &hitout.;
   set hits;
   length blast_hit_flag $32.;
   blast_hit_flag="&flag.";
run;

%mend;

%extractHits(match_flags,flag_len_query_eq_ref,true_hit,match_flags2);
* 14655 true hits;

%extractHits(match_flags2,flag_len_query_eq_ref_mm5,true_hit_mm5,match_flags3);
* 620 true hits with 5 MM max;

%extractHits(match_flags3,flag_len_match_eq_ref,complete_hit,match_flags4);
* 21445 complete hits;

%extractHits(match_flags4,flag_len_match_eq_ref_mm5,complete_hit_mm5,match_flags5);
* 163 complete hits with 5 MM max;

%extractHits(match_flags5,flag_len_match_5bp_diff_ref,hit_diff_5bp,match_flags6);
* 8684 hits where match is within 5bp of reference length;

%extractHits(match_flags6,flag_len_match_5bp_diff_ref_mm5,hit_diff_5bp_mm5,match_flags7);
* 155 hits where match is within 5bp of reference length and 5 MM max;

%extractHits(match_flags7,flag_len_match_10bp_diff_ref,hit_diff_10bp,match_flags8);
* 2959 hits where match is within 10bp of reference length;

%extractHits(match_flags8,flag_len_match_10bp_diff_ref_mm5,hit_diff_10bp_mm5,match_flags9);
* 100 hits where match is within 10bp of reference length and 5 MM max;

%extractHits(match_flags9,flag_len_match_90prc_ref,match_90perc,match_flags10);
* 25358 hits;

%extractHits(match_flags10,flag_len_match_90prc_ref_mm5,match_90perc_mm5,no_hits);
* 738 hits;
* in total 418779 "non-hits";


/* Stack hits and save */

data event.sim_xscript2refseq_blast_best;
   set true_hit true_hit_mm5 complete_hit complete_hit_mm5
       hit_diff_5bp hit_diff_5bp_mm5 hit_diff_10bp hit_diff_10bp_mm5
       match_90perc match_90perc_mm5;
run;

proc sort data=event.sim_xscript2refseq_blast_best ;
   by simulation test blast_hit_flag;
proc freq data=event.sim_xscript2refseq_blast_best noprint;
   by simulation test;
   tables blast_hit_flag / out=blast_summary;
run;

proc print data=blast_summary;
run;

proc freq data=event.sim_xscript2refseq_blast_best ;
   tables simulation*test;
run;

/*
      Table of simulation by test

  simulation     test

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |test1   |test2   |  Total
  ---------+--------+--------+
  sim1     |  12686 |  12434 |  25120
           |  16.94 |  16.61 |  33.55
           |  50.50 |  49.50 |
           |  33.70 |  33.40 |
  ---------+--------+--------+
  sim2     |  12508 |  12447 |  24955
           |  16.70 |  16.62 |  33.33
           |  50.12 |  49.88 |
           |  33.23 |  33.43 |
  ---------+--------+--------+
  sim3     |  12451 |  12351 |  24802
           |  16.63 |  16.50 |  33.12
           |  50.20 |  49.80 |
           |  33.07 |  33.17 |
  ---------+--------+--------+
  Total       37645    37232    74877
              50.28    49.72   100.00
*/









