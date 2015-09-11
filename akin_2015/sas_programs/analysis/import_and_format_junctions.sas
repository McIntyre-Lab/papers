/* Because the junctions datasets are large I will be importing and analyzing the data by chromosome, then stacking later */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';

/* import junctions */

    data WORK.JUNCTION_COVERAGE    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "/scratch/lfs/sugrue/splicing_by_chrom/all_counts_&sysparm..csv"
 delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $19. ;
       informat fusion_id $2475. ;
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
       format sample_id $19. ;
       format fusion_id $2475. ;
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

*called events as fusions. I will need to fix this up later;


/* summarize then annotate */

*drop unneeded variables - rpkm, mean, std, cv, reads_in_region, mapped_reads;

data junction_coverage;
   set junction_coverage;
   drop rpkm mean std cv read_length reads_in_region mapped_reads;
run;


/* need to add subject_id to counts file */
/* then sum mapped reads, region depth, reads in region, and apn by subject and gene */
/* can drop rpkm, mean, std, cv */

proc sort data=mysas.subject2sample;
by sample_id;
run;


proc sort data=junction_coverage;
by sample_id;
run;

data jnc_counts_w_key oops1 oops2;
   merge mysas.subject2sample (in=in1) junction_coverage (in=in2);
   by sample_id;
   if in1 and in2 then output jnc_counts_w_key;
   else if in1 then output oops1;
   else output oops2;
run;


proc datasets noprint;
   delete junction_coverage;
run; quit;


/* split counts_w_key region_depth, reads_in_region, apn */


proc sort data=jnc_counts_w_key;
   by subject fusion_id;
run;

proc means data=jnc_counts_w_key noprint;
   var region_depth;
   by subject fusion_id;
   output out=jnc_region_depth sum=;
run;

proc means data=jnc_counts_w_key noprint;
   var apn;
   by subject fusion_id;
   output out=jnc_apn sum=;
run;


/* merge all together */

data jnc_subject_info;
   set jnc_counts_w_key;
   keep subject fusion_id region_length;
run;

proc sort data=jnc_subject_info nodup;
  by subject fusion_id region_length;
run;

proc sort data=jnc_region_depth;
by subject fusion_id;
run;

proc sort data=jnc_apn;
by subject fusion_id;
run;

data jnc_subject_depth oops1 oops2;
  merge jnc_subject_info (in=in1) jnc_region_depth (in=in2);
  by subject fusion_id;
  if in1 and in2 then output jnc_subject_depth;
  else if in1 then output oops1;
  else output oops2;
 drop _TYPE_ _FREQ_;
run;

data mysas.jnc_subject_apn_&sysparm. oops1 oops2;
  merge jnc_subject_depth (in=in1) jnc_apn (in=in2);
  by subject fusion_id;
  if in1 and in2 then output mysas.jnc_subject_apn_&sysparm.;
  else if in1 then output oops1;
  else output oops2;
 drop _TYPE_ _FREQ_;
run;

/* add in treatment */

proc sort data=mysas.sample_key;
by subject;
run;

proc sort data=mysas.jnc_subject_apn_&sysparm.;
by subject;
run;

data subject_counts_w_key oops1 oops2;
   merge mysas.jnc_subject_apn_&sysparm. (in=in1) mysas.sample_key (in=in2);
   by subject;
   if in1 and in2 then output subject_counts_w_key;
   else if in1 then output oops1;
   else output oops2;
run;


/* add in log-apn and make permenant */

data mysas.jnc_counts_w_key_&sysparm.;
   set subject_counts_w_key;
   apn=region_depth/region_length;
   log_apn=log(apn+1);
   if apn>0 then flag_junction_on=1;
   else flag_junction_on=0;
run;





