/* 
 * REVISIONS: 12/21/2011 
 *            - Changed some folder structure
 *            - Had to change to full import statment because sample_id was not
 *              importing properly            
 */

libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';


 data ALL_COVERAGE_COUNTS    ;
 %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
 infile '!MCLAB/arbeitman/arbeitman_fru_network/pipeline_output/coverage_on_fusions/all_coverage_counts.csv' 
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



 data truncated    ;
 %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
 infile '!MCLAB/arbeitman/arbeitman_fru_network/pipeline_output/coverage_on_fusions/coverage_on_fusions_2011-07-05_3_TTAGGC.csv'
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=1 ;
    informat fusion_id $9. ;
    informat sample_id $20. ;
    informat tmapped_reads best32. ;
    informat treads_in_exon best32. ;
    informat tcoverage_in_exon best32. ;
    informat exon_length best32. ;
    informat tapn best32. ;
    informat trpkm best32. ;
    format fusion_id $9. ;
    format sample_id $20. ;
    format tmapped_reads best12. ;
    format treads_in_exon best12. ;
    format tcoverage_in_exon best12. ;
    format exon_length best12. ;
    format tapn best12. ;
    format trpkm best12. ;
 input
             fusion_id $
             sample_id
             tmapped_reads
             treads_in_exon
             tcoverage_in_exon
             exon_length
             tapn
             trpkm
 ;
 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;


 data miss    ;
 %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
 infile '!MCLAB/arbeitman/arbeitman_fru_network/pipeline_output/coverage_on_fusions/coverage_on_fusions_2011-07-05_3_AH_Male_FruM_A_3_missing.csv'
    delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=1 ;
    informat fusion_id $9. ;
    informat msample_id $20. ;
    informat mmapped_reads best32. ;
    informat mreads_in_exon best32. ;
    informat mcoverage_in_exon best32. ;
    informat mexon_length best32. ;
    informat mapn best32. ;
    informat mrpkm best32. ;
    format fusion_id $9. ;
    format msample_id $20. ;
    format mmapped_reads best12. ;
    format mreads_in_exon best12. ;
    format mcoverage_in_exon best12. ;
    format mexon_length best12. ;
    format mapn best12. ;
    format mrpkm best12. ;
 input
             fusion_id $
             msample_id
             mmapped_reads
             mreads_in_exon
             mcoverage_in_exon
             mexon_length
             mapn
             mrpkm
 ;
 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;



proc sort data=truncated;
    by fusion_id;
    run;

proc sort data=miss;
    by fusion_id;
    run;

data merged oops;
    merge truncated (in=in1) miss (in=in2);
    by fusion_id;
    if in1 and in2 then output merged;
    else output oops;
    run;

data recalc;
    retain fusion_id sample_id mapped_reads reads_in_exon coverage_in_exon exon_length apn rpkm;
    set merged;
    mapped_reads = tmapped_reads + mmapped_reads;
    coverage_in_exon = tcoverage_in_exon + mcoverage_in_exon;
    reads_in_exon = coverage_in_exon / 70;
    apn = coverage_in_exon / exon_length;
    rpkm = 1000000000 * reads_in_exon / (mapped_reads * exon_length);
    keep fusion_id sample_id mapped_reads reads_in_exon coverage_in_exon exon_length apn rpkm;
    run;
    

data clean;
set all_coverage_counts ;
if sample_id eq '2011-07-05_3_TTAGGC' then delete;
run;

data fru.all_coverage_counts_miss;
set clean recalc;
run; * 3617460 obs;




