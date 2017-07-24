ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';


/* Going to try UQ normalizing the TPM counts to see if that helps */

*macro-tized so I can use the same code;


* Do this on only the set of samples I've been using;

data subject_list;
   set eventloc.subjects_for_rsem_test;
run;

data design;
   set con.design_by_subject_new;
   keep subject_id cell_type name library;
run;

proc sort data=subject_list;
   by subject_id;
proc sort data=design;
   by subject_id;
run;

data design2;
   merge design (in=in1) subject_list (in=in2);
   by subject_id;
   if in1 and in2;
run;

%macro norm(datain,dataout);

data counts;
  set eventloc.&datain.;
run;

proc sort data=counts;
   by library;
proc sort data=design2;
   by library;
run;

data counts_w_key;
  merge counts (in=in1) design2 (in=in2);
  by library;
  if in1 and in2;
run;

data counts_w_key2;
   set counts_w_key;
   if tpm > 0;
run;

proc univariate data=counts_w_key2 noprint;
   by library;
   var tpm;
   output out=quartiles_tpm0 Q3=Q3;
run;

proc means data=quartiles_tpm0 noprint;
    var Q3;
    output out=q3_q3 q3=q3;
run;

%local upQuart;

data _null_;
   set q3_q3;
   call symputx('upQuart',Q3);
run;


proc sort data=counts_w_key2;
   by library;
proc sort data=quartiles_tpm0;
   by library;
run;

data counts_w_key_q3;
  merge counts_w_key2 (in=in1) quartiles_tpm0 (in=in2);
  by library;
  if in1 and in2;
run;

data &dataout.;
   set counts_w_key_q3;
   log_tpm=log(tpm+1);
   q3_tpm=(tpm/q3) * &upQuart. ;
   q3_ff=&upQuart. / q3 ;
   log_q3_tpm=log(q3_tpm + 1);
run;

%mend;

%norm(hg19_rsem_all_xscripts,hg19_rsem_all_xs_q3norm);
%norm(hg19_rsem_75perc_apn5_xscripts,hg19_rsem_75perc_apn5_xs_q3norm);

data hg19_rsem_all_xs_q3norm2;
  set hg19_rsem_all_xs_q3norm;
  keep library transcript_id name cell_type subject_id log_q3_tpm;
run;

data hg19_rsem_75perc_apn5_xs_q3norm2;
  set hg19_rsem_75perc_apn5_xs_q3norm;
  keep library transcript_id name cell_type subject_id log_q3_tpm;
run;

/* Make permenant */

data eventloc.hg19_rsem_all_xs_q3norm;
  set hg19_rsem_all_xs_q3norm;
run;

data eventloc.hg19_rsem_75perc_apn5_xs_q3norm;
  set hg19_rsem_75perc_apn5_xs_q3norm;
run;



proc export data=hg19_rsem_all_xs_q3norm2 outfile="/mnt/store/event_sandbox/tpm_by_subject_cell_all.csv"
  dbms=csv replace;
run;

proc export data=hg19_rsem_75perc_apn5_xs_q3norm2 outfile="/mnt/store/event_sandbox/tpm_by_subject_cell_filtered.csv"
  dbms=csv replace;
run;

