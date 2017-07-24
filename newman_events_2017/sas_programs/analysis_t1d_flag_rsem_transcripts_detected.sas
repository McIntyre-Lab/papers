ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For each transcript list, flag transcripts that are "on"
   e.g. if TPM>0 then transcript is detected
   then it must be "detected" in at least 50% of subjects or "not detected" in any */

data subject_list;
   set eventloc.subjects_for_rsem_test;
run;

data design;
   set con.design_by_subject_new;
   keep subject_id cell_type library;
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


%macro flagXSon(datain,dataout,minTPM);

data counts;
  set eventloc.&datain.;
  if tpm >= &minTPM. then flag_tpm_gt0=1;
  else flag_tpm_gt0=0;
run;

proc sort data=counts;
    by library;
proc sort data=design2;
    by library;
run;

data counts2;
  merge counts (in=in1) design2 (in=in2);
  by library;
  if in1 and in2;
run;

proc sort data=counts2;
   by cell_type transcript_id;
run;

proc means data=counts2 noprint;
   by cell_type transcript_id;
   var flag_tpm_gt0;
   output out=xs_on_by_trt mean=perc_trt_on;
run;

data cd4 cd8 cd19;
   set xs_on_by_trt;
   if cell_type='CD4' then output cd4;
   if cell_type='CD8' then output cd8;
   if cell_type='CD19' then output cd19;
   keep transcript_id perc_trt_on;
run;

data cd4_2;
  set cd4;
  if perc_trt_on ge 0.5 then flag_cd4_on=1;
  else if perc_trt_on=0 then flag_cd4_on=0;
  else flag_cd4_on=.;
  keep transcript_id flag_cd4_on;
run;

data cd8_2;
  set cd8;
  if perc_trt_on ge 0.5 then flag_cd8_on=1;
  else if perc_trt_on=0 then flag_cd8_on=0;
  else flag_cd8_on=.;
  keep transcript_id flag_cd8_on;
run;

data cd19_2;
  set cd19;
  if perc_trt_on ge 0.5 then flag_cd19_on=1;
  else if perc_trt_on=0 then flag_cd19_on=0;
  else flag_cd19_on=.;
  keep transcript_id flag_cd19_on;
run;

proc sort data=cd4_2;
   by transcript_id;
proc sort data=cd8_2;
   by transcript_id;
proc sort data=cd19_2;
   by transcript_id;
run;

data &dataout.;
   merge cd4_2 cd8_2 cd19_2;
   by transcript_id;
run;
%mend;

%flagXSon(hg19_rsem_all_xscripts,flag_xscript_tpm_ge1_all,1);
%flagXSon(hg19_rsem_75perc_apn5_xscripts,flag_xscript_tpm_ge1_filt,1);
%flagXSon(hg19_rsem_all_xscripts,flag_xscript_tpm_ge5_all,5);
%flagXSon(hg19_rsem_75perc_apn5_xscripts,flag_xscript_tpm_ge5_filt,5);

/* Make permenant */

data eventloc.flag_xscript_tpm_ge1_all;
  set flag_xscript_tpm_ge1_all;
run;

data eventloc.flag_xscript_tpm_ge1_filt;
   set flag_xscript_tpm_ge1_filt;
run;

data eventloc.flag_xscript_tpm_ge5_all;
  set flag_xscript_tpm_ge5_all;
run;

data eventloc.flag_xscript_tpm_ge5_filt;
   set flag_xscript_tpm_ge5_filt;
run;


* check detection overlap;

proc freq noprint data=flag_xscript_tpm_ge1_all;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=xs_count;
proc print data=xs_count;
run;


/*
  flag_     flag_     flag_
 cd4_on    cd8_on    cd19_on     COUNT
   .         .         .        38558
   .         .         0        25250
   .         .         1         3842
   .         0         .         4704
   .         0         0        10841
   .         0         1          254
   .         1         .         1752
   .         1         0          520
   .         1         1         1332
   0         .         .         5967
   0         .         0        46617
   0         .         1          604
   0         0         .        11504
   0         0         0       386848
   0         0         1         1342
   0         1         .           29
   0         1         0          122
   0         1         1           27
   1         .         .         2651
   1         .         0          505
   1         .         1         1190
   1         0         .            5
   1         0         0           38
   1         0         1            4
   1         1         .         6899
   1         1         0         1262
   1         1         1        43420


*/

proc freq noprint data=flag_xscript_tpm_ge1_filt;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=xs_count;
proc print data=xs_count;
run;


/*
  flag_     flag_     flag_
 cd4_on    cd8_on    cd19_on    COUNT

    .         .         .        5480
    .         .         0        1283
    .         .         1         991
    .         0         .         548
    .         0         0         707
    .         0         1          46
    .         1         .         412
    .         1         0          80
    .         1         1         386
    0         .         .         616
    0         .         0        1003
    0         .         1         142
    0         0         .         982
    0         0         0        5721
    0         0         1         268
    0         1         .           2
    0         1         0          17
    0         1         1           8
    1         .         .         451
    1         .         0          65
    1         .         1         315
    1         0         .           1
    1         0         0           8
    1         0         1           2
    1         1         .        1459
    1         1         0         263
    1         1         1       17408

*/




proc freq noprint data=flag_xscript_tpm_ge5_all;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=xs_count;
proc print data=xs_count;
run;


/*
  flag_     flag_     flag_
 cd4_on    cd8_on    cd19_on    COUNT
    .         .         .         6820
    .         .         0         6112
    .         .         1         1473
    .         0         .         1048
    .         0         0         2850
    .         0         1          112
    .         1         .          400
    .         1         0          242
    .         1         1          466
    0         .         .         1304
    0         .         0        30429
    0         .         1          297
    0         0         .         4168
    0         0         0       521794
    0         0         1          792
    0         1         .           14
    0         1         0           77
    0         1         1           25
    1         .         .          503
    1         .         0          292
    1         .         1          344
    1         0         .            6
    1         0         0           22
    1         0         1            3
    1         1         .         1744
    1         1         0         762
    1         1         1       13988


*/


proc freq noprint data=flag_xscript_tpm_ge5_filt;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=xs_count;
proc print data=xs_count;
run;

/*
  flag_     flag_     flag_
 cd4_on    cd8_on    cd19_on    COUNT
    .         .         .        3005
    .         .         0        1689
    .         .         1         687
    .         0         .         312
    .         0         0         629
    .         0         1          49
    .         1         .         241
    .         1         0         108
    .         1         1         225
    0         .         .         441
    0         .         0        1932
    0         .         1         123
    0         0         .        1020
    0         0         0       17708
    0         0         1         331
    0         1         .           3
    0         1         0          39
    0         1         1           9
    1         .         .         259
    1         .         0         108
    1         .         1         180
    1         0         0          13
    1         0         1           1
    1         1         .        1010
    1         1         0         397
    1         1         1        8145
*/


