
libname relapse "/home/ammorse/staph_relapse/sasdata";
libname lmm "/home/ammorse/staph_relapse/sasdata/data_analysis";


/*

use sra_sampleID from staph.AF_kitchen_sink_locked to 'rename' fq files for submission to SRA

output: 	staph.design_sra_submission  (MCLAB/staph/staph_DNA_v2/design_files/mscramms/design_sra_submission.csv")
			columns:  isolate  database_submissionID  read  fq  set


*/

proc contents data = relapse.sampleID_st_cc_ref_kitchen_sink ; run;

data meta ;
set relapse.sampleID_st_cc_ref_kitchen_sink ;
where flag_failed_library = 0 ;
organism = "Staphylococcus aureus" ;
strain = "missing";
isolate = "missing";
host = "Homo sapiens";
host_disease = "Staphylococcus aureus bacteremia";
collected_by = "SAB Group Prospective Cohort Study at Duke University" ;
geo_loc_name = "USA" ;
isolation_source = "blood culture" ;
lat_lon = "missing" ;
keep sra_sampleID organism strain isolate host collected_by geo_loc_name isolation_source sampleID patientnumber host_disease lat_lon;
run ;

/* use BC date for collection date */
proc import datafile = "/home/ammorse/staph_relapse/Recur_SAB_isolates_Data_Codebook_2018_05_18_Lauren3-17-20_isolatedCheckedOut_sheet.csv"
out = dating 
dbms = csv replace ;
run;

data bc_date ;
set dating ;
format collection_date DATE9. ;
rename study_number = sampleID ;
collection_date = BC_date ;
keep study_number patient_number collection_date ;
run;

proc sort data = meta ;
by sampleID ;
proc sort data = bc_date ;
by sampleID ;
run;

data meta_w_bcDate ;
merge meta (in=in1) bc_date (in=in2) ;
by sampleID ;
if in1 ;
run;

data check ;
set meta_w_bcDate  ;
where patientNumber ne patient_number ;
run;  /* 0 obs */

data relapse.design_sra_attrib ;
retain sra_sampleID organism strain isolate collected_by collection_date geo_loc_name host host_disease isolation_source lat_lon host_subject_id ;
set meta_w_bcDate ;
host_subject_id = patientnumber ;
drop patientnumber sampleID patient_number ;
label sra_sampleID = "sample_name" ;
run ;

proc freq data =  relapse.design_sra_attrib ;
tables  sra_sampleID / out = cnts ;
run ;
data check ;
set cnts ;
where count ne 1 ;
run;

proc export data = relapse.design_sra_attrib 
outfile = "/home/ammorse/staph_relapse/sra/design_sra_attrib.tsv"
label 
dbms = tab replace ;
run;


/* create file with 'meta data' -- includes fq info, etc
Sample name (db submission ID)
Library ID (same as sampleName) 
title, WGS of S. aureus isolate XXXXX
Library strategy (WGS) 
Library source (GENOMIC), 
Library selection (RANDOM) 
Library layout paired
Platform ILLUMINA
Instrument model Illumina HiSeq2500
Design description brief M&M
Filetype fastq
Filename 
Filename1

*/

data ID_list;
set relapse.sampleID_st_cc_ref_kitchen_sink ;
keep sampleID sra_sampleID ;
run;

data fq_list ;
set relapse.fq_list_w_sampleNum_sampleID ;
where flag_failed_library = 0 ;
keep fqName read sampleID;
run;

proc sort data = ID_list ;
by sampleID ;
proc sort data = fq_list ;
by sampleID ;
run ;

data df_fq missing;
merge ID_list (in=in1) fq_list (in=in2) ;
by sampleID ;
if in2  then output df_fq ;
else output missing ;
run ;  /* 1 in missing is the failed library */


/* create file with 'meta data' -- includes fq info, etc
Sample name (db submission ID)
Library ID (same as sampleName) 
title, WGS of S. aureus isolate XXXXX
Library strategy (WGS) 
Library source (GENOMIC), 
Library selection (RANDOM) 
Library layout paired
Platform ILLUMINA
Instrument model Illumina HiSeq2500
Design description brief M&M
Filetype fastq
Filename 
Filename1

*/

data df_fq2 ;
set df_fq ;
/* add needed submission columns */
title = ("WGS of S. aureus sample "||sra_sampleID) ;
library_strategy = "WGS";
library_source = "GENOMIC";
library_selection = "RANDOM";
library_layout = "paired" ;
platform = "ILLUMINA";
instrument_model = "Illumina NovaSeq 6000";
design_description = "DNA was purified from each sample followed by whole-genome library preparation and sequencing. Barcoded libraries were pooled and 2 x 150bp reads were generated on an Illumina NovaSeq 6000." ;
filetype = "fastq";  
sampleID ;
run ;

proc sort data = df_fq2 ;
by sra_sampleID ;
run ;

proc transpose data = df_fq2 out = wide ;
by sra_sampleID ;
var fqName ;
run ;

data reads ;
set wide ;
rename col1 = filename ;
rename col2 = filename2 ;
drop _name_ ;
run ;

data df_fq3 ;
set df_fq2;
drop fqName read ;
run ;

proc sort data = df_fq3 nodups ;
by _all_ ;
run;

proc sort data = df_fq3 ;
by sra_sampleID ;
proc sort data = reads ;
by sra_sampleID ;
run ;

data df_sra_close fix;
merge df_fq3 (in=in1) reads (in=in2) ;
by sra_sampleID ;
if in1 and in2 then output df_sra_close;
else output fix ;
run ;  /* 0 in fix */

data df_sra_ready ;
retain sra_sampleID library_ID title library_strategy library_source library_selection library_layout platform instrument_model design_description filetype filename filename2 ;
set  df_sra_close ;
library_ID = sra_sampleID ;
label sra_sampleID = "sample_name" ;
run ;

data relapse.design_sra_meta ;
set df_sra_ready ;
run;


proc sort data = relapse.design_sra_meta  nodups dupout= dups;
by filename filename2 ;
run;  /* as expected, no dups */

proc sort data = relapse.design_sra_meta  ;
by sra_sampleID  ;
run;

proc export data = relapse.design_sra_meta 
outfile = "/home/ammorse/staph_relapse/sra/design_sra_meta.tsv"
label
dbms = tab replace ;
run;


/* design for copying fq files for packing - tar files in a list (omit failed library) */

data read1 ;
set relapse.design_sra_meta  ;
keep sra_sampleID filename ;
run ;

data read2 ;
set relapse.design_sra_meta  ;
keep sra_sampleID filename2 ;
rename filename2 = filename ;
run ;

data relapse.design_sra_pack_fq ;
set read1 read2 ;
drop sra_sampleID ;
run ;

proc sort data = relapse.design_sra_pack_fq ;
by fileName  ;
run;

proc export data = relapse.design_sra_pack_fq 
outfile = "/home/ammorse/staph_relapse/sra/design_sra_pack_fq.csv"
dbms = csv replace ;
run;


