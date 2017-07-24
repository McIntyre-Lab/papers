ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Import PacBio BLAST results -- I am going stack together the fragment and junction BLAST results and process together */


    data WORK.PACBIO_FRAG    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/analysis_output/blast_output/blast_fragments_to_pacbio.tsv'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat feature_id $173. ;
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
        format feature_id $173. ;
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
                feature_id $
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


    data WORK.PACBIO_JUNC    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'!MCLAB/event_analysis/analysis_output/blast_output/blast_events_to_pacbio_min_length_50.tsv'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat feature_id $173. ;
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
       format feature_id $173. ;
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
                feature_id $
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


data pacbio_blast_results;
   length feature_type $8.;
   set pacbio_frag (in=in1) pacbio_junc (in=in2);
   if in1 then feature_type="fragment";
   if in2 then feature_type="splicing";
run;

/* Remove hits with gaps or mismatches */

data pacbio_blast_results2;
  set pacbio_blast_results;
  if mismatch > 0 then delete;
  if gapopen > 0 then delete;
run;

/* Split PacBio ID */

data pacbio_blast_results3;
   length pacbio_gene_id $15.;
   length pacbio_xscript_id $15.;
   length pacbio_status $10.;
   length pacbio_type $50.;
   set pacbio_blast_results2;
   pacbio_xscript_id=scan(pacbio_id,1,"|");
   pacbio_status=scan(pacbio_id,2,"|");
   pacbio_type=scan(pacbio_id,3,"|");
   pacbio_gene_id=catx(".","PB",scan(pacbio_xscript_id,2,"."));
   drop pacbio_id;
   rename pacbio_xscript_id=pacbio_id;
run;

/* Make permenant */

data event.pacbio_blast_results;
   set pacbio_blast_results3;
run;




