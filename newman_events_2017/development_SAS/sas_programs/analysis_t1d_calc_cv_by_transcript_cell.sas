ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For each transcript list, calculate the variance in TPM per transcript per cell type */

data design;
   set con.design_by_subject_new;
   /* Drop these samples -- low coverage */
   if name =  '2009-PC-0221' then delete; *sample 75 cd8;
   if name =  '2009-PC-0144' then delete; *sample 48 cd4;
   if name =  '2009-PC-0236' then delete; *sample 80;
   if name =  '2009-PC-0237' then delete; *sample 80;
   if name =  '2009-PC-0235' then delete; *sample 80; 
   keep subject_id cell_type library;
run;

%macro calcVar(datain,dataout);

data counts;
  set eventloc.&datain.;
  log_tpm=log(tpm+1);
run;

proc sort data=counts;
   by library;
proc sort data=design;
   by library;
run;

data counts_w_key;
  merge counts (in=in1) design (in=in2);
   by library;
  if in1 and in2;
run;

proc sort data=counts_w_key;
  by cell_type transcript_id;
proc means data=counts_w_key noprint;
  by cell_type transcript_id;
  var log_tpm;
  output out=var_by_xs_cell mean=tpm_mean var=tpm_var cv=tpm_cv;
run;

proc sort data=var_by_xs_cell;
   by transcript_id cell_type;
proc transpose data=var_by_xs_cell out=var_sbys;
   by transcript_id;
   var tpm_var;
   id cell_type;
run;

proc transpose data=var_by_xs_cell out=mean_sbys;
   by transcript_id;
   var tpm_mean;
   id cell_type;
run;

proc transpose data=var_by_xs_cell out=cv_sbys;
   by transcript_id;
   var tpm_cv;
   id cell_type;
run;

data var_sbys2;
   set var_sbys;
   keep transcript_id cd4 cd8 cd19;
   rename cd4=var_cd4 cd8=var_cd8 cd19=var_cd19;
run;

data mean_sbys2;
   set mean_sbys;
   keep transcript_id cd4 cd8 cd19;
   rename cd4=mean_cd4 cd8=mean_cd8 cd19=mean_cd19;
run;

data cv_sbys2;
   set cv_sbys;
   keep transcript_id cd4 cd8 cd19;
   rename cd4=cv_cd4 cd8=cv_cd8 cd19=cv_cd19;
run;

proc sort data=mean_sbys2;
   by transcript_id;
proc sort data=var_sbys2;
   by transcript_id;
proc sort data=cv_sbys2;
   by transcript_id;
run;


data &dataout. ;
  merge mean_sbys2 (in=in1) var_sbys2 (in=in2) cv_sbys2 (in=in2);
   by transcript_id;
  if in1 and in2;
run;

%mend;

%calcVar(hg19_rsem_all_xscripts,hg19_cv_all_xs);
%calcVar(hg19_rsem_75perc_apn5_xscripts,hg19_cv_filt_xs);

/* Make permenant */

data eventloc.hg19_cv_all_xs;
   set hg19_cv_all_xs;
run;

data eventloc.hg19_cv_filtered_xs;
   set hg19_cv_filt_xs;
run;

