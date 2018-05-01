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
  set event.ireckon2pacbio_blast;
run;

proc sort data=blast_hits;
   by sample query_id transcript_id;
proc freq data=blast_hits noprint;
   by sample query_id;
   tables transcript_id / out=ref_hits_per_query;
run;

data fragmented;
  set ref_hits_per_query;
  where count > 1;
  keep sample query_id transcript_id;
run;

proc sort data=blast_hits;
   by sample query_id transcript_id;
proc sort data=fragmented;
   by sample query_id transcript_id;
run;

data blast_hits_no_frag;
  merge blast_hits (in=in1) fragmented (in=in2);
  by sample query_id transcript_id;
  if in1 and in2 then delete;
  else if in1 then output;
run; *92954 blast hits in total, reduced to 6934;

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
  keep sample query_id;
run;

proc sort data=others;
  by sample query_id;
proc sort data=query2drop nodup;
  by sample query_id;
run;

data &otherout.;
  merge others (in=in1) query2drop (in=in2);
  by sample query_id;
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
* 5 true hits;

%extractHits(match_flags2,flag_len_query_eq_ref_mm5,true_hit_mm5,match_flags3);
* 1 true hit with 5 MM max;

%extractHits(match_flags3,flag_len_match_eq_ref,complete_hit,match_flags4);
* 2203 complete hits;

%extractHits(match_flags4,flag_len_match_eq_ref_mm5,complete_hit_mm5,match_flags5);
* 19 complete hits with 5 MM max;

%extractHits(match_flags5,flag_len_match_5bp_diff_ref,hit_diff_5bp,match_flags6);
* 65 hits where match is within 5bp of reference length;

%extractHits(match_flags6,flag_len_match_5bp_diff_ref_mm5,hit_diff_5bp_mm5,match_flags7);
* 7 hits where match is within 5bp of reference length and 5 MM max;

%extractHits(match_flags7,flag_len_match_10bp_diff_ref,hit_diff_10bp,match_flags8);
* 22 hits where match is within 10bp of reference length;

%extractHits(match_flags8,flag_len_match_10bp_diff_ref_mm5,hit_diff_10bp_mm5,match_flags9);
* 1 hit where match is within 10bp of reference length and 5 MM max;

%extractHits(match_flags9,flag_len_match_90prc_ref,match_90perc,match_flags10);
* 567 hits;

%extractHits(match_flags10,flag_len_match_90prc_ref_mm5,match_90perc_mm5,no_hits);
* 16 hits;
* in total 4028 "non-hits";


/* Stack hits and save */

data event.ireckon2pacbio_blast_best;
   set true_hit true_hit_mm5 complete_hit complete_hit_mm5
       hit_diff_5bp hit_diff_5bp_mm5 hit_diff_10bp hit_diff_10bp_mm5
       match_90perc match_90perc_mm5;
run;

proc sort data=event.ireckon2pacbio_blast_best ;
   by sample blast_hit_flag;
proc freq data=event.ireckon2pacbio_blast_best noprint;
   by sample;
   tables blast_hit_flag / out=blast_summary;
run;

proc print data=blast_summary;
run;
/*
sample    blast_hit_flag                      COUNT    PERCENT

 NSC1     flag_len_match_10bp_diff_ref           12     0.6893
 NSC1     flag_len_match_10bp_diff_ref_mm5        1     0.0574
 NSC1     flag_len_match_5bp_diff_ref            43     2.4698
 NSC1     flag_len_match_5bp_diff_ref_mm5         2     0.1149
 NSC1     flag_len_match_90prc_ref              311    17.8633
 NSC1     flag_len_match_90prc_ref_mm5            8     0.4595
 NSC1     flag_len_match_eq_ref                1352    77.6565
 NSC1     flag_len_match_eq_ref_mm5              11     0.6318
 NSC1     flag_len_query_eq_ref                   1     0.0574
 NSC2     flag_len_match_10bp_diff_ref           10     0.8584
 NSC2     flag_len_match_5bp_diff_ref            22     1.8884
 NSC2     flag_len_match_5bp_diff_ref_mm5         5     0.4292
 NSC2     flag_len_match_90prc_ref              256    21.9742
 NSC2     flag_len_match_90prc_ref_mm5            8     0.6867
 NSC2     flag_len_match_eq_ref                 851    73.0472
 NSC2     flag_len_match_eq_ref_mm5               8     0.6867
 NSC2     flag_len_query_eq_ref                   4     0.3433
 NSC2     flag_len_query_eq_ref_mm5               1     0.0858

*/

proc freq data=event.ireckon2pacbio_blast_best ;
   tables sample;
run;

/*

                                     Cumulative    Cumulative
  sample    Frequency     Percent     Frequency      Percent
  -----------------------------------------------------------
  NSC1          1741       59.91          1741        59.91
  NSC2          1165       40.09          2906       100.00


*/









