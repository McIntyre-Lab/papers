/* Import coverage counts. Need all counts in one enormous file. Need to add a header */
/* Make stacked dataset with key*/

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data2';

/*First import the coverage counts after combining all 5024 files*/
/*Will make a sas dataset on my local rather than running everything through mclab*/
/* 1:11 to import file into sas, maybe there is a way to subset the data to make this go quicker */

  data mysas.all_counts_by_gene_&sysparm.    ;
      %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
      infile "/scratch/lfs/patcon/jnewman/splicing_by_chrom_cell/all_counts_&sysparm..csv" delimiter = ',' MISSOVER
 DSD lrecl=32767 firstobs=2;
         informat sample_id $33. ;
         informat event_id $2475. ;
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
         format sample_id $33. ;
         format event_id $2475. ;
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





proc sort data=mysas.all_counts_by_gene_&sysparm.;
  by sample_id;
  run;

data design_file_&sysparm.;
    set mysas.design_file_se;
run;


proc sort data=design_file_&sysparm.;
  by sample_id;
  run;


data mysas.all_counts_w_key_&sysparm.;
  merge mysas.all_counts_by_gene_&sysparm. (in=in1) design_file_&sysparm. (in=in2);
  by sample_id;
  *logapn = log(apn +1);
  if in1;
run;

/*proc export data=mysas.all_counts_w_key
    outfile='/scratch/lfs/patcon/jnewman/sas_analysis/all_counts_splicing_w_key.csv'
    label dbms=csv replace;
    run;
*/



*Sort by transcript_id;
proc sort data=mysas.all_counts_w_key_&sysparm.;
  by Name event_id;
  run;

*Sum technical reps for each sample;
proc means data = mysas.all_counts_w_key_&sysparm. noprint;
  *class Name;
  by Name event_id;
  var region_depth;
  output out=mysas.mapped_reads_summed_byevent_&sysparm. sum=depth;
run;

proc datasets noprint;
   delete mysas.all_counts_w_key_&sysparm.;
run;
quit;

* Recalculate APN and make dataset permanent;
/*data mysas.counts_by_event;
   set mysas.mapped_reads_summed_byevent;
   if name ne " " ;
   apn = depth / region_length;
   run;
*/
* calculate total mapped reads by subject;
/*
data mysas.mapped_reads ;
  set mysas.all_counts_w_key;
  keep subject_id mapped_reads;
  run;

proc sort data=mysas.mapped_reads nodup;
  by subject_id mapped_reads;
  run;

proc means data=mysas.mapped_reads;
  by subject_id;
  var mapped_reads;
  output out=mysas.total_mapped_reads_sum sum=total_mapped_reads;
  run;
*/

*export for check;
/*proc export data=mysas.mapped_reads_summed_byevent
    outfile='/scratch/lfs/patcon/jnewman/sas_analysis/mapped_reads_summed.csv'
    label dbms=csv replace;
    run;

proc export data=mysas.total_mapped_reads_sum
    outfile='/scratch/lfs/patcon/jnewman/sas_analysis/total_mapped_reads.csv'
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
