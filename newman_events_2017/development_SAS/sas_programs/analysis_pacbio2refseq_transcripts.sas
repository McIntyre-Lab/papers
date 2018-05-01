ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Match PacBio transcripts to RefSeq transcripts */

/* Import splice-match transcripts */

     data WORK.SPLICE_MATCH    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/references/pacbio_isoforms_list_for_import.txt'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat pacbio_id $10. ; informat source $6.;    informat feature_type $10. ;
        informat start best32. ;  informat end best32.;   informat frame best32. ;
        informat strand $1. ;     informat score best32.; informat transcript_id $18. ;
        informat match_type $23.; informat note $28. ;
        format pacbio_id $10. ;   format source $6. ;     format feature_type $10. ;
        format start best12. ;    format end best12. ;    format frame best12. ;
        format strand $1. ;       format score best12. ;  format transcript_id $18. ;
        format match_type $23. ;  format note $28. ;
        input pacbio_id $ source $ feature_type $ start end
              frame strand $ score transcript_id $ match_type $ note $
              ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

/* Import BLAST results */

proc import datafile="!MCLAB/event_analysis/analysis_output/blast_output/known_pacbio2refseq_best_hits.tsv"
      out=blast_hits dbms=tab replace;
     guessingrows=11000; getnames=no;
run;

    data WORK.BLAST_HITS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/event_analysis/analysis_output/blast_output/known_pacbio2refseq_best_hits.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat query_id $40. ; informat transcript_id $12. ; informat perc_identity best32. ;
       informat length best32. ; informat mismatch best32. ; informat gapopen best32. ;
       informat query_start best32. ; informat query_stop best32. ; informat ref_start best32. ;
       informat ref_stop best32. ; informat evalue best32. ; informat bitscore best32. ;
       format query_id $40. ;  format transcript_id $12. ;  format perc_identity best12. ;
       format length best12. ; format mismatch best12. ; format gapopen best12. ;
       format query_start best12. ; format query_stop best12. ;  format ref_start best12. ;
       format ref_stop best12. ;  format evalue best12. ; format bitscore best12. ;
       input query_id $ transcript_id $ perc_identity length mismatch gapopen
             query_start query_stop ref_start ref_stop evalue bitscore ;
  if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
  run;


/* Get the list of Refseq matches from the splice-match set */

data splice_match2;
   set splice_match;
   where transcript_id ^? "ENSMUST";
   length splice_match_id $18.;
   splice_match_id=scan(transcript_id,1,".");
   keep pacbio_id splice_match_id;
run;

/* Process BLAST hits */

proc sort data=blast_hits;
   by query_id;
run;


data blast_hits2;
   length pacbio_id $10.;
   length match_type $30.;
   set blast_hits;
   by query_id;
   pacbio_id=scan(query_id,1,"|");
   match_type=scan(query_id,3,"|");
   if first.query_id then output;
   keep pacbio_id match_type transcript_id;
   rename transcript_id=blast_hit_id;
run;

proc sort data=blast_hits2;
   by pacbio_id;
proc sort data=splice_match2;
   by pacbio_id;
run;

data splice_match2blast_hit;
  merge splice_match2 (in=in1) blast_hits2 (in=in2);
  by pacbio_id;
run;

/* Check concordance between splice match and BLAST hit
   If different, take the splice match. If no splice match, take the BLAST hit */

data flag_diff;
   set splice_match2blast_hit;
   if splice_match_id="" then flag_refseq_diff=-1;
   else if blast_hit_id="" then flag_refseq_diff=-2;
   else if splice_match_id=blast_hit_id then flag_refseq_diff=0;
   else flag_refseq_diff=1;
run;


proc freq data=flag_diff;
   tables flag_refseq_diff;
run;


/*

flag_refseq_diff    Frequency     Percent
-------------------------------------------
              -2           1        0.01
              -1         931        9.64
               0        7837       81.19
               1         884        9.16

1 instance where there is no BLAST hit
931 instances of no splice match (Ensembl transcripts?)
7837 instances where splice match and BLAST hit are the same
884 instances where BLAST hit is different to splice hit
*/

data set_refseq_transcript;
   set flag_diff;
   length transcript_id $30.;
   if flag_refseq_diff=-1 then transcript_id=blast_hit_id;
   else transcript_id=splice_match_id;
   if transcript_id="" then flag_no_id=1;
   else flag_no_id=0;
run;

proc freq data=set_refseq_transcript;
   tables flag_no_id;
run;

/*
 flag_no_id    Frequency     Percent
 ------------------------------------
          0        9653      100.00


all have a refseq ID!
*/

/* Make permenant */

data event.pacbio2refseq_id;
  set set_refseq_transcript;
  keep pacbio_id transcript_id match_type;
run;



