ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For each set of RSEM estimates, I want to count the number of subjects in each cell type */

data subj_all;
   set eventloc.hg19_rsem_all_xscripts;
   keep library;
run;

data subj_filt;
   set eventloc.hg19_rsem_75perc_apn5_xscripts;
   keep library;
run;

data lib2cell;
  set con.design_by_subject_new;
  keep library subject_id name cell_type;
run;

proc sort data=subj_all nodup;
  by library;
proc sort data=subj_filt nodup;
  by library;
proc sort data=lib2cell nodup;
   by library subject_id name cell_type;
run;

data subj_all_w_cell;
  merge subj_all (in=in1) lib2cell (in=in2);
  by library;
  if in1 and in2;
run;

data subj_filt_w_cell;
  merge subj_filt (in=in1) lib2cell (in=in2);
  by library;
  if in1 and in2;
run;

proc freq data=subj_filt_w_cell;
  tables cell_type;
run;

/*
    cell_type    Frequency
    -------------------------
    CD19               30
    CD4                30
    CD8                30
*/

proc freq data=subj_all_w_cell;
  tables cell_type;
run;

/*

 cell_type    Frequency
 --------------------------
 CD19               31
 CD4                32
 CD8                29

Okay, so both ~30
*/

*How many are common to both sets?;

proc sort data=subj_filt_w_cell;
   by name;
proc sort data=subj_all_w_cell;
   by name;
run;

data subject_check;
  merge subj_filt_w_cell (in=in1) subj_all_w_cell (in=in2);
  by name;
  if in1 then flag_in_filt=1; else flag_in_filt=0;
  if in2 then flag_in_all=1; else flag_in_all=0;
run;

proc freq data=subject_check noprint;
   tables cell_type*flag_in_filt*flag_in_all / out=check;
proc print data=check;
run;

/*
   cell_    flag_in_    flag_in_
   type       filt         all      COUNT    PERCENT

   CD19         0           1          4         4
   CD19         1           0          3         3
   CD19         1           1         27        27
   CD4          0           1          3         3
   CD4          1           0          1         1
   CD4          1           1         29        29
   CD8          0           1          3         3
   CD8          1           0          4         4
   CD8          1           1         26        26



Okay, of the ones that are in both, how many subjects have all three in both lists? */

data cd4 cd8 cd19;
   set subject_check;
   if flag_in_filt=1 and flag_in_all=1 then do;
     if cell_type="CD4" then output cd4;
     if cell_type="CD8" then output cd8;
     if cell_type="CD19" then output cd19;
   end;
   keep subject_id;
run;

proc sort data=cd4;
  by subject_id;
proc sort data=cd8;
  by subject_id;
proc sort data=cd19;
  by subject_id;
run;

data check2;
  merge cd4 (in=in1) cd8 (in=in2) cd19 (in=in3);
  by subject_id;
  if in1 then flag_cd4=1; else flag_cd4=0;
  if in2 then flag_cd8=1; else flag_cd8=0;
  if in3 then flag_cd19=1; else flag_cd19=0;
run;

proc freq data=check2 noprint;
  tables flag_cd4*flag_cd8*flag_cd19 / out=check_count;
proc print data=check_count;
run;

/*
                         flag_
 flag_cd4    flag_cd8     cd19    COUNT

     0           1         1         1
     1           0         0         1
     1           0         1         3
     1           1         0         2
     1           1         1        23

okay subset this list of 23 subjects -- will use this for testing
*/

data eventloc.subjects_for_rsem_test;
  set check2;
  where flag_cd4=1 and flag_cd8=1 and flag_cd19=1;
  keep subject_id;
run;

