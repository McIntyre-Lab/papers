/* import libraries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* import sugrue gene coverage data for ANOVAs */

proc import datafile='/home/jrbnewman/McLab/sugrue/original_data/coverage_counts_fusions_strand.csv'
   out=counts
   dbms=csv replace;
   getnames=yes;
   guessingrows=370000;
run;

/* need to add subject_id to counts file */
/* then sum mapped reads, region depth, reads in region, and apn by subject and gene */
/* can drop rpkm, mean, std, cv */

proc sort data=sugrue.subject2sample;
by sample_id;
run;

proc sort data=counts;
by sample_id;
run;

data counts_w_key oops1 oops2;
   merge sugrue.subject2sample (in=in1) counts (in=in2);
   by sample_id;
   if in1 and in2 then output counts_w_key;
   else if in1 then output oops1;
   else output oops2;
run;

/* split counts_w_key region_depth, reads_in_region, apn */

proc sort data=counts_w_key;
   by subject fusion_id;
run;

proc means data=counts_w_key noprint;
   var region_depth;
   by subject fusion_id;
   output out=region_depth sum=;
run;

proc means data=counts_w_key noprint;
   var reads_in_region;
   by subject fusion_id;
   output out=reads_in_region sum=;
run;

proc means data=counts_w_key noprint;
   var apn;
   by subject fusion_id;
   output out=apn sum=;
run;

/* merge all together */

data subject_info;
   set counts_w_key;
   keep subject fusion_id read_length region_length;
run;

proc sort data=subject_info nodup;
  by subject fusion_id read_length region_length;
run;

proc sort data=region_depth;
by subject fusion_id;
run;

proc sort data=reads_in_region;
by subject fusion_id;
run;

proc sort data=apn;
by subject fusion_id;
run;

data subject_depth oops1 oops2;
  merge subject_info (in=in1) region_depth (in=in2);
  by subject fusion_id;
  if in1 and in2 then output subject_depth;
  else if in1 then output oops1;
  else output oops2;
 drop _TYPE_ _FREQ_;
run;

data subject_reads oops1 oops2;
  merge subject_depth (in=in1) reads_in_region (in=in2);
  by subject fusion_id;
  if in1 and in2 then output subject_reads;
  else if in1 then output oops1;
  else output oops2;
 drop _TYPE_ _FREQ_;
run;

data subject_apn oops1 oops2;
  merge subject_reads (in=in1) apn (in=in2);
  by subject fusion_id;
  if in1 and in2 then output subject_apn;
  else if in1 then output oops1;
  else output oops2;
 drop _TYPE_ _FREQ_;
run;

/* add in treatment */

proc sort data=sugrue.sample_key;
by subject;
run;

proc sort data=subject_apn;
by subject;
run;

data subject_counts_w_key oops1 oops2;
   merge subject_apn (in=in1) sugrue.sample_key (in=in2);
   by subject;
   if in1 and in2 then output subject_counts_w_key;
   else if in1 then output oops1;
   else output oops2;
run;


/* add in log-apn and make permenant */

data sugrue.counts_w_key;
   set subject_counts_w_key;
   apn=region_depth/region_length;
   log_apn=log(apn+1);
   if apn>0 then flag_fusion_on=1;
   else flag_fusion_on=0;
run;


