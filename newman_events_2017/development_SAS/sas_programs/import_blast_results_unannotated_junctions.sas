ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Import BLAST results for unannotated junctions */

    data WORK.UNANNOT_BLAST    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/event_analysis/analysis_output/blast_output/blast_unannotated_junctions_to_pacbio_basic.tsv' delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat event_id $146. ;
       informat pacbio_id $40. ;
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
       format event_id $146. ;
       format pacbio_id $40. ;
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
                event_id $
                pacbio_id $
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

*2527 alignments;

/* Check : how many alignments have <80% identity? */

data flag_ident;
   set unannot_blast;
   if perc_identity < 90 then flag_perc_identity_lt_90=1;
   else flag_perc_identity_lt_90=0;
run;

proc freq data=flag_ident;
   tables flag_perc_identity_lt_90;
run; *0 < 90%!;


/* Check : how many alignments have a length <50 bp? */

data flag_length;
   set unannot_blast;
   if length < 50 then flag_length_lt_50=1;
   else flag_length_lt_50=0;
run;

proc freq data=flag_length;
   tables flag_length_lt_50;
run;

*0 alignments less than 50 bp;

/* Split Pacbio ID into transcript_id, known or novel, match type */

data unannot_blast_parse_pb;
  length pb_transcript_id $12.;
  length pb_status $8.;
  length pb_type $30.;
  set unannot_blast;
  pb_transcript_id=scan(pacbio_id,1,"|");
  pb_status=scan(pacbio_id,2,"|");
  pb_type=scan(pacbio_id,3,"|");
  query_length=abs(query_stop-query_start)+1;
  ref_length=abs(ref_stop-ref_start)+1;
run;

data flag_lengths;
   set unannot_blast_parse_pb;
   if query_length = ref_length then flag_query_length_eq_ref=1;
   else flag_query_length_eq_ref=0;

   if query_length = length then flag_query_length_eq_len=1;
   else flag_query_length_eq_len=0;

   if ref_length = length then flag_ref_length_eq_len=1;
   else flag_ref_length_eq_len=0;
run;

proc freq data=flag_lengths noprint;
  tables flag_query_length_eq_ref*flag_query_length_eq_len*flag_ref_length_eq_len / out=length_check;
run;

proc print data=length_check;
run;

/*
 flag_query_    flag_query_    flag_ref_
   length_        length_       length_
    eq_ref         eq_len        eq_len     COUNT

     0              0             0            1     0.0396
     0              0             1           61     2.4139
     0              1             0          123     4.8674
     1              0             0            5     0.1979
     1              1             1         2337    92.4812


Most of these have the query and ref lengths the same
*/

/* Flag if perc identity > 90%, and length >= 45

   We have event sizes of 80bp
   For hits of 40 bp or less, these could be matches to only half of the junction (ie, only the donor, or only the acceptor)
   If we use 41 bp as the minimum, then this means that the alignment should span the exon-exon junction
   and hence will be the "best hits". If we also filter out alignments <90%, then we should also get the "best hits".
   The overlap between these two variables will be the set of "confirmed junctions".

   If any junction does not have a "good hit", then take the top hit of each, as we only want to see if the junction
   matches to any novel transcript */

data flag_best_hits;
   set unannot_blast_parse_pb;
   if perc_identity ge 90 then flag_perc_identity_ge_90=1;
   else flag_perc_identity_ge_90=0;
   if length ge 41 then flag_length_ge_41bp=1;
   else flag_length_ge_41bp=0;
   if mismatch > 0 then flag_mismatch=1; else flag_mismatch=0;
   if gapopen > 0 then flag_gapopen=1; else flag_gapopen=0;
run;


proc freq data=flag_best_hits noprint;
   tables flag_perc_identity_ge_90*flag_length_ge_41bp*flag_mismatch*flag_gapopen / out=blast_check;
run;

proc print data=blast_check; run;

/* Want to keep only those hits where the perc identity is >90%, the length is at least 45bp
  and there are no gaps (mismatches are okay)
 flag_perc_
  identity_    flag_length_      flag_      flag_
    ge_90         ge_41bp      mismatch    gapopen    COUNT    PERCENT	KEEP?

      1              1             0          0        2274    89.9881	Y
      1              1             0          1         160     6.3316	N
      1              1             1          0          63     2.4931	Y?
      1              1             1          1          30     1.1872	N


2274 "good" hits, plus 63 possible extra

*/

data unannot_best_hits;
   set flag_best_hits;
   if flag_perc_identity_ge_90=1 and flag_length_ge_41bp=1 and flag_gapopen=0;
run;

data count_uniq_events;
   set unannot_best_hits;
   keep event_id;
run;

proc sort data=count_uniq_events nodup;
  by event_id;
run; *1470 unique hits;

/* Make permenant */

data event.unannot_junc_best_pacbio_hits;
   set unannot_best_hits;
run;

data event.unannot_junc_pacbio_hits;
   set flag_best_hits;
run;


proc freq data=event.unannot_junc_best_pacbio_hits;
   table pb_type;
run;


/*

  pb_type                 Frequency     Percent
  ----------------------------------------------
  antisense                         1        0.04
  full-splice_match              1308       55.97
  fusion                            5        0.21
  genic                             5        0.21
  genic_intron                      2        0.09
  incomplete-splice_match         120        5.13
  intergenic                        1        0.04
  novel_in_catalog                634       27.13
  novel_not_in_catalog            261       11.17


Lots matching to "full-splice_match" transcripts...

*/

/* Check lengths of hits against lengths of junction.
   Subset only the hits where the BLAST length is the same as the junction length, and there are no gaps
   Then count:
   novel/known * flag_junction_annotated

   Also compare against catalog! (ie, merge PB junctions to catalog junctions on coordinate,
   how many hit catalog junctions? What is the crossover with the PB type here? */
  



