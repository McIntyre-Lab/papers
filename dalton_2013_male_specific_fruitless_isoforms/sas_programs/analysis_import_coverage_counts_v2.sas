/* 
 * REVISIONS: 12/21/2011 
 *            - Changed some folder structure
 *            - Had to change to full import statment because sample_id was not
 *              importing properly            
 */

libname fru '!MCLAB/Fru_network/sasdata';


 data FRU.ALL_COVERAGE_COUNTS    ;
 %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
 infile '!MCLAB/Fru_network/pipeline_output/coverage_on_fusions/all_coverage_counts.csv' 
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
    informat fusion_id $9. ;
    informat sample_id $20. ;
    informat mapped_reads best32. ;
    informat reads_in_exon best32. ;
    informat coverage_in_exon best32. ;
    informat exon_length best32. ;
    informat apn best32. ;
    informat rpkm best32. ;
    format fusion_id $9. ;
    format sample_id $20. ;
    format mapped_reads best12. ;
    format reads_in_exon best12. ;
    format coverage_in_exon best12. ;
    format exon_length best12. ;
    format apn best12. ;
    format rpkm best12. ;
 input
             fusion_id $
             sample_id
             mapped_reads
             reads_in_exon
             coverage_in_exon
             exon_length
             apn
             rpkm
 ;
 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;
