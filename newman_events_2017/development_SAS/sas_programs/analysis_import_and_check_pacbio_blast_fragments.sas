ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Import fragment PacBio BLAST hits and process */

    data WORK.PACBIO_FRAG    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/analysis_output/blast_output/blast_fragments_to_pacbio.tsv'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat fragment_id $13. ;
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
        format fragment_id $13. ;
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
                fragment_id $
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
proc sort data=pacbio_frag;
   by fragment_id;
run;

data pacbio_frag_w_flags;
  merge pacbio_frag (in=in1) frag_on frag_info frag_short;
  by fragment_id;
  if in1;
run;


/* Parse PacBio ID, so that I can map these back to RefSeq transcripts */

data pacbio_frag_w_flags2;
   length pacbio_id2 $15.;
   set pacbio_frag_w_flags;
   pacbio_id2=scan(pacbio_id,1,"|");
   drop pacbio_id;
   rename pacbio_id2=pacbio_id;
run;


/* Keeping hits where:
   (1) percent identity is 100%
   (2) alignmnet length = fragment length
   (3) no gaps or mismatches
   (4) fragment length is at least 12bp

What I can do to check on 3' variability is flag fragments from the 3'-most exon of the transcript,
   and see if all match, or if only a subset match

i.e.

Transcript: ----------[===============================]
Fragments:  ----------[======|========|=========|=====]
Hits:       ----------[======]
Detection   ----------[======|========]

Could be used to indicate 3' variation

*/

data  hits_kept hits_dropped;
   set pacbio_frag_w_flags2;
   if flag_fragment_lt_min_bp=1 then delete;
   if perc_identity = 100 and length = fragment_length
   and gapopen = 0 and mismatch = 0 then output hits_kept;
   else output hits_dropped;
run;


/* PacBio 2 RefSeq */

data pb2refseq;
   set event.pacbio2refseq_id;
   keep pacbio_id transcript_id;
run;

proc sort data=pb2refseq nodup;
   by pacbio_id;
proc sort data=hits_kept;
   by pacbio_id;
run;

data hits_pb2refseq;
   merge hits_kept (in=in1) pb2refseq (in=in2);
   by pacbio_id;
   if in2 then flag_has_refseq_id=1; else flag_has_refseq_id=0;
   if in1 then output;
run;

/* Check if hits are going to PacBio-counterparts of their original ref seq transcripts */

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
proc sort data=hits_pb2refseq;
   by fragment_id transcript_id;
run;

data hits_w_orig_xs hits_w_other_xs hits_w_no_xs;
    merge hits_pb2refseq (in=in1) frag2xs (in=in2);
   by fragment_id transcript_id;
   if in1 and in2 then output hits_w_orig_xs;
   else if in1 and transcript_id="" then output hits_w_no_xs;
   else if in1 and transcript_id ne "" then output hits_w_other_xs;
run;

/* Make output permenant for now. Will need to think about these later */

data event.frag_pacbio_blast_hits_known;
   set hits_w_orig_xs;
run;

data event.frag_pacbio_blast_hits_other;
   set hits_w_other_xs;
run;

data event.frag_pacbio_blast_hits_no_refseq;
   set hits_w_no_xs;
run;




