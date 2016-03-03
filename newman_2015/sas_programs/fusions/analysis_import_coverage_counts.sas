/* Import coverage counts. Need all counts in one enormous file. Need to add a header */
/* Make stacked dataset with key*/

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data';

/*First import the coverage counts after combining all 5024 files*/
/*Will make a sas dataset on my local rather than running everything through mclab*/
/* 1:11 to import file into sas, maybe there is a way to subset the data to make this go quicker */

  data mysas.all_counts_by_gene    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile '/scratch/lfs/patcon/jnewman/all_counts_by_fusion.csv' delimiter = ',' MISSOVER
 DSD lrecl=32767 firstobs=2;
         informat sample_id $31. ;
         informat fusion_id $24. ;
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
         format sample_id $31. ;
         format fusion_id $24. ;
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





proc sort data=mysas.all_counts_by_gene;
  by sample_id;
  run;


proc sort data=mysas.design_file;
  by sample_id;
  run;


data mysas.all_counts_w_key;
  merge mysas.all_counts_by_gene (in=in1) mysas.design_file (in=in2);
  by sample_id;
  *logapn = log(apn +1);
  if in1;
run;
/*
proc export data=mysas.all_counts_w_key
    outfile='/scratch/lfs/patcon/jnewman/sas_analysis/all_counts_byfusion_w_key.csv'
    label dbms=csv replace;
    run;
*/




/* Going to make a test set to run through all of the scripts	 *
 * Subset 500 transcripts and merge with design file 		 */
/* Figure out how many obs will give you 500 transcripts

proc sort data=mysas.tech_reps_summed;
  by transcript_id;
  run;

data test;
  set mysas.tech_reps_summed;
  obs=_n_;
  if obs le 1500;
  run;
*/
