/* First import the coverage counts. Then create a stacked dataset of all the counts with the design file
*/


libname "/home/fnew/mclab/arbeitman/arbeitman_ribotag/sas_data";

data work.all_coverage_counts ;
 %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile
 '/home/fnew/mclab/arbeitman/arbeitman_ribotag/pipeline_output/coverage_count_uniq_2/all_cvg_cnts_uniq.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
         informat sample_id $28. ;
         informat fusion_id $9. ;
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
         format sample_id $28. ;
         format fusion_id $9. ;
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

data ribo.all_coverage_counts2;
set all_coverage_counts;
run;


/* Now merge with design file to get stacked dataset */

proc sort data=ribo.all_coverage_counts2;
by sample_id;
run;

proc sort data=ribo.design_file;
by sample_id;
run;

data ribo.all_coverage_counts_with_key2;
merge ribo.all_coverage_counts2 (in=in1) ribo.design_file (in=in2);
by sample_id ;
logrpkm = log(rpkm +1) ;
if in1;
run ;
*1389982 observations;

