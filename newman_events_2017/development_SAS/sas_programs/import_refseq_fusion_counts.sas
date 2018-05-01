libname event '!MCLAB/event_analysis/sas_data';

/* Import fusion counts */

    data WORK.COUNTS_BY_FUSION    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/alignment_output/mm10_refseq_fusion_counts.csv' delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $4. ;
       informat fusion_id $13. ;
       informat mapped_reads best32. ;
       informat read_length best32. ;
       informat region_length best32. ;
       informat region_depth best32. ;
       informat reads_in_region best32. ;
       informat apn best32. ;
       informat rpkm best32. ;
       informat mean best32. ;
       informat std best32. ;
       informat cv best32. ;
       format sample_id $4. ;
       format fusion_id $13. ;
       format mapped_reads best12. ;
       format read_length best12. ;
       format region_length best12. ;
       format region_depth best12. ;
       format reads_in_region best12. ;
       format apn best12. ;
       format rpkm best12. ;
       format mean best12. ;
       format std best12. ;
       format cv best12. ;
    input
                sample_id $
                fusion_id $
                mapped_reads
                read_length
                region_length
                region_depth
                reads_in_region
                apn
                rpkm
                mean
                std
                cv
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Make permenant */

data event.mm10_refseq_fusion_counts;
  set counts_by_fusion;
run;

