ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Import BLAST results for splicing events to PacBio transcripts and RefSeq transcripts */

/* PacBio BLAST */

    data WORK.PACBIO_BLAST    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/event_analysis/analysis_output/blast_output/blast_events_to_pacbio_min_length_50.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat event_id $173. ;
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
       format event_id $173. ;
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

/* RefSeq BLAST */

    data WORK.REFSEQ_BLAST    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/event_analysis/analysis_output/blast_output/blast_events_to_refseq_min_length_50.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat event_id $173. ;
       informat refseq_id $12. ;
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
       format event_id $173. ;
       format refseq_id $12. ;
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
                refseq_id $
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


/* Make permenant */

data event.dtct_event_pacbio_blast_results;
   set pacbio_blast;
run;

data event.dtct_event_refseq_blast_results;
   set refseq_blast;
run;

