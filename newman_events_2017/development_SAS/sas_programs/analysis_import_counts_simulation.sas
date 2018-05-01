/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Import simulated counts and format */

/* Import splicing counts */

    data WORK.SPLICING_COUNTS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/alignment_output/simulations/mm10_refseq_splicing_counts_simulation.csv' delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $10. ;
       informat event_id $335. ;
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
       format sample_id $10. ;
       format event_id $335. ;
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
               event_id $
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

/* Fusion counts */

    data WORK.COUNTS_BY_FUSION    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/alignment_output/simulations/mm10_refseq_fusion_counts_simulation.csv' delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $10. ;
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
       format sample_id $10. ;
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


/* Fragment counts */
    data WORK.COUNTS_BY_FRAGMENT    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '!MCLAB/event_analysis/alignment_output/simulations/mm10_refseq_fragment_counts_simulation.csv' delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $10. ;
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
       format sample_id $10. ;
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
data event.mm10_refseq_splicing_counts_sim;
 set SPLICING_COUNTS;
 keep sample_id event_id region_depth apn;
run;

data event.mm10_refseq_fusion_counts_sim;
  set counts_by_fusion;
  keep sample_id fusion_id apn;
run;

data event.mm10_refseq_fragment_counts_sim;
  set counts_by_fragment;
  keep sample_id fusion_id apn;
  rename fusion_id=fragment_id;
run;
