ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Import fragment RefSeq BLAST hits and process */

    data WORK.REFSEQ_FRAG    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/analysis_output/blast_output/blast_fragments_to_refseq.tsv'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat fragment_id $13. ;
        informat transcript_id $12. ;
        informat perc_identity best32. ;
        informat length best32. ;
        informat mismatch best32. ;
        informat gapopen best32. ;
        informat query_start best32. ;
        informat query_stop best32. ;
        informat ref_start best32. ;
        informat ref_stop best32. ;
        informat evalue best32. ;
        informat bitscore best32. ;
        format fragment_id $13. ;
        format transcript_id $12. ;
       format perc_identity best12. ;
       format length best12. ;
       format mismatch best12. ;
       format gapopen best12. ;
       format query_start best12. ;
       format query_stop best12. ;
       format ref_start best12. ;
       format ref_stop best12. ;
       format evalue best12. ;
       format bitscore best12. ;


     input
                fragment_id $
                transcript_id $
                 perc_identity
                length
                mismatch
                gapopen
                query_start
                query_stop
                ref_start
                ref_stop
                evalue
                bitscore
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Get length and detection flags */

data frag_on;
   set event.fragments_on_apn_gt0;
run;

data frag_short;
   set event.flagged_fragment_length;
run;

data frag_info;
   set mm10.mm10_exon_fragment_flagged;
   fragment_length=fragment_end-fragment_start;
   keep fragment_id fragment_length transcript_id;
   rename transcript_id=transcript_id_cat;
run;

proc sort data=frag_on;
   by fragment_id;
proc sort data=frag_short;
   by fragment_id;
proc sort data=frag_info;
   by fragment_id;
proc sort data=refseq_frag;
   by fragment_id;
run;

data refseq_frag_w_flags;
  merge refseq_frag (in=in1) frag_on frag_info frag_short;
  by fragment_id;
  if in1;
run;


/* Because we are checking if fragments are being assigned to the transcripts they have come from,
   and no fragments are unannotated, alignment are kept if:
   (1) percent identity is 100%
   (2) alignmnet length = fragment length
   (3) no gaps or mismatches
   (4) fragment length is at least 12bp */

data  hits_kept hits_dropped;
   set refseq_frag_w_flags;
   if flag_fragment_lt_min_bp=1 then delete;
   if perc_identity = 100 and length = fragment_length
   and gapopen = 0 and mismatch = 0 then output hits_kept;
   else output hits_dropped;
run;

proc means data=hits_kept noprint;
   var evalue;
   output out=evalue_dist mean=mean stddev=sd max=max min=min
       q1=q1 median=median q3=q3;
run;

proc print data=evalue_dist;run;

/* 
mean=3.706385E-21
sd=4.439342E-20
min= 0
q1= 4E-103
median=8E-64
q3=2E-42
max= 7E-19
*/


data hits_kept2;
  set hits_kept;
  keep fragment_id;
run;

data hits_dropped2;
  set hits_dropped;
  keep fragment_id;
run;

proc sort data=hits_kept2 nodup;
   by fragment_id;
proc sort data=hits_dropped2 nodup;
   by fragment_id;
run;

data hits_check;
  merge hits_kept2 (in=in1) hits_dropped2 (in=in2);
  by fragment_id;
  if in1 then flag_frag_hit=1; else flag_frag_hit=0;
  if in2 then flag_frag_no_hit=1; else flag_frag_no_hit=0;
run;0

proc freq data=hits_check;
  tables flag_frag_hit*flag_frag_no_hit;
run;

/*

  flag_frag_hit     flag_frag_no_hit

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |   3247 |   3247
           |   0.00 |   2.06 |   2.06
           |   0.00 | 100.00 |
           |   0.00 |  25.78 |
  ---------+--------+--------+
         1 | 144836 |   9347 | 154183
           |  92.00 |   5.94 |  97.94
           |  93.94 |   6.06 |
           | 100.00 |  74.22 |
  ---------+--------+--------+
  Total      144836    12594   157430
              92.00     8.00   100.00

Got some fragments with no hit...
*/

/* For hits, check that they map to the same set of transcripts from which they are derived */


data frag2xs;
   set mm10.mm10_exon_fragment_flagged;
   length orig_transcript_id $20.;
   do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      orig_transcript_id=scan(transcript_id,i,"|");
      output;
      end;
   keep fragment_id orig_transcript_id;
   rename orig_transcript_id=transcript_id;
run;

proc sort data=frag2xs nodup;
   by fragment_id transcript_id;
proc sort data=hits_kept;
   by fragment_id transcript_id;
run;

data hits_w_orig_xs hits_w_other_xs;
    merge hits_kept (in=in1) frag2xs (in=in2);
   by fragment_id transcript_id;
   if in1 and in2 then output hits_w_orig_xs;
   else if in1 then output hits_w_other_xs;
run;
*16394 fragments with alternative hits;

/* how many are below the median evalue of original hits? */


proc means data=hits_w_orig_xs noprint;
   var evalue;
   output out=evalue_dist mean=mean stddev=sd max=max min=min
       q1=q1 median=median q3=q3;
run;

proc print data=evalue_dist;run;

/* 
mean=3.706385E-21
sd=4.439342E-20
min= 0
q1= 4E-104
median=6E-65
q3=2E-42
max= 7E-19
*/


data hits_other_flag_med;
   set hits_w_other_xs;
   if evalue < 6E-65 then flag_hit_lt_med=1;
   else flag_hit_lt_med=0;
run;

proc freq data=hits_other_flag_med;
  tables flag_hit_lt_med;
run;

/*
                                               Cumulative    Cumulative
   flag_hit_lt_med    Frequency     Percent     Frequency      Percent
   --------------------------------------------------------------------
                 0       12307       75.07         12307        75.07
                 1        4087       24.93         16394       100.00

*/


/* Make output permenant for now. Will need to think about these later */

data event.frag_refseq_blast_hits_known;
   set hits_w_orig_xs;
run;

data event.frag_refseq_blast_hits_other;
   set hits_w_other_xs;
run;




