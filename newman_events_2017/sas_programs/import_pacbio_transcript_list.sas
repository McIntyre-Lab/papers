ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';

/* Import list of PacBio transcripts */

    data WORK.PB_XS_LIST    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/references/pacbio_transcript_list.txt' delimiter='09'x
MISSOVER DSD lrecl=32767 ;
       informat pacbio_id $10. ;
       informat pb_status $5. ;
       informat pb_type $23. ;
       format pacbio_id $10. ;
       format pb_status $5. ;
       format pb_type $23. ;
    input
                pacbio_id $
                pb_status $
                pb_type $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Add gene ID */

data pb_xs_w_gene;
   length pacbio_gene_id $10.;
   set pb_xs_list;
   pacbio_gene_id=catt("PB.",scan(pacbio_id,2,"."));
run;

/* Make permenant */

data event.pacbio_transcripts;
  set pb_xs_w_gene;
run;

