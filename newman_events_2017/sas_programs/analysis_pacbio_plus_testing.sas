
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';

/* Waiting on PacBio-plus RSEM estimates, so I am going to "simulate" these from the existing PB results, and just add in the missing transcripts for the other */

data list1 list2;
  set event.xscripts_w_unique_by_bin;
  if perc_features_dtct >= 0.5 then output list1;
  if perc_features_dtct >= 0.75 then output list2;
   keep transcript_id;
run;

data list3 list4;
  set event.bin_xscripts_by_dtct_apn5;
  if perc_features_dtct >= 0.5 then output list3;
  if perc_features_dtct >= 0.75 then output list4;
   keep transcript_id;
run;

proc sort data=list1;
   by transcript_id;
proc sort data=list2;
   by transcript_id;
proc sort data=list3;
   by transcript_id;
proc sort data=list4;
   by transcript_id;
run;

data xs_list;
  merge list1 (in=in1) list2 (in=in2) list3 (in=in3) list4 (in=in4);
  by transcript_id;
  if in1 then xs_50perc_apn0=1; else xs_50perc_apn0=0;
  if in2 then xs_75perc_apn0=1; else xs_75perc_apn0=0;
  if in3 then xs_50perc_apn5=1; else xs_50perc_apn5=0;
  if in4 then xs_75perc_apn5=1; else xs_75perc_apn5=0;
run;


/* Get estimates from previous RSEM rounds */

data est1;
  set event.rsem_events_exp_50perc;
  keep transcript_id tpm_nsc1 tpm_nsc2;
  rename tpm_nsc1=tpm_nsc1_50p_a0 tpm_nsc2=tpm_nsc2_50p_a0;
run;

data est2;
  set event.rsem_events_exp_75perc;
  keep transcript_id tpm_nsc1 tpm_nsc2;
  rename tpm_nsc1=tpm_nsc1_75p_a0 tpm_nsc2=tpm_nsc2_75p_a0;
run;

data est3;
  set event.rsem_events_exp_50perc_apn5;
  keep transcript_id tpm_nsc1 tpm_nsc2;
  rename tpm_nsc1=tpm_nsc1_50p_a5 tpm_nsc2=tpm_nsc2_50p_a5;
run;

data est4;
  set event.rsem_events_exp_75perc_apn5;
  keep transcript_id tpm_nsc1 tpm_nsc2;
  rename tpm_nsc1=tpm_nsc1_75p_a5 tpm_nsc2=tpm_nsc2_75p_a5;
run;


proc sort data=est1;
  by transcript_id;
proc sort data=est2;
  by transcript_id;
proc sort data=est3;
  by transcript_id;
proc sort data=est4;
  by transcript_id;
proc sort data=xs_list;
  by transcript_id;
run;

data xs_w_est;
  merge xs_list (in=in1) est1 est2 est3 est4;
  by transcript_id;
  if in1 ;
run;

/* Add in PacBio results */

data xs_w_pb;
   set event.pacbio2refseq_id;
   keep transcript_id pacbio_id;
run;

data pb_est;
   set event.rsem_pacbio_all;
   keep transcript_id tpm_nsc1 tpm_nsc2;
   rename transcript_id=pacbio_id tpm_nsc1=tpm_nsc1_pb tpm_nsc2=tpm_nsc2_pb;
run;

proc sort data=xs_w_pb;
   by pacbio_id;
proc sort data=pb_est;
   by pacbio_id;
run;

data xs2pb_est;
  merge xs_w_pb (in=in1) pb_est (in=in2);
   by pacbio_id;
  if in2;
run;

proc sort data=xs2pb_est;
   by transcript_id;
proc sort data=xs_w_est nodup;
   by transcript_id;
run;

data xs_list_w_pb;
  merge xs2pb_est (in=in1) xs_w_est (in=in2);
  by transcript_id;
  if in1 then do;
    xs_50perc_apn0=1; xs_75perc_apn0=1; xs_50perc_apn5=1; xs_75perc_apn5=1;
  end;
run;

/* Subset lists for plots */

data pb_plus_50p_apn0;
   set xs_list_w_pb;
   where xs_50perc_apn0=1;
   if pacbio_id ne "" then do;
          tpm_nsc1=tpm_nsc1_pb;
          tpm_nsc2=tpm_nsc2_pb;
          end;
   else do;
          tpm_nsc1=tpm_nsc1_50p_a0;
          tpm_nsc2=tpm_nsc2_50p_a0;
          end;
   keep pacbio_id transcript_id tpm_nsc1 tpm_nsc2;
run;


data pb_plus_75p_apn0;
   set xs_list_w_pb;
   where xs_75perc_apn0=1;
   if pacbio_id ne "" then do;
          tpm_nsc1=tpm_nsc1_pb;
          tpm_nsc2=tpm_nsc2_pb;
          end;
   else do;
          tpm_nsc1=tpm_nsc1_75p_a0;
          tpm_nsc2=tpm_nsc2_75p_a0;
          end;
   keep pacbio_id transcript_id tpm_nsc1 tpm_nsc2;
run;

data pb_plus_50p_apn5;
   set xs_list_w_pb;
   where xs_50perc_apn5=1;
   if pacbio_id ne "" then do;
          tpm_nsc1=tpm_nsc1_pb;
          tpm_nsc2=tpm_nsc2_pb;
          end;
   else do;
          tpm_nsc1=tpm_nsc1_50p_a5;
          tpm_nsc2=tpm_nsc2_50p_a5;
          end;
   keep pacbio_id transcript_id tpm_nsc1 tpm_nsc2;
run;


data pb_plus_75p_apn5;
   set xs_list_w_pb;
   where xs_75perc_apn5=1;
   if pacbio_id ne "" then do;
          tpm_nsc1=tpm_nsc1_pb;
          tpm_nsc2=tpm_nsc2_pb;
          end;
   else do;
          tpm_nsc1=tpm_nsc1_75p_a5;
          tpm_nsc2=tpm_nsc2_75p_a5;
          end;
   keep pacbio_id transcript_id tpm_nsc1 tpm_nsc2;
run;

/* Bin and calc rep concordance */

%macro concord(xslist);

data tpm_data;
   set &xslist.;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id pacbio_id mean_tpm tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
run;

proc means data=tpm_data noprint;
  where mean_tpm > 0;
  var mean_tpm;
  output out=tpm_distrib min=min q1=q1 q3=q3 max=max;
run;

%local tpmMin tpmQ1 tpmQ3 tpmMax;

data _null_;
   set tpm_distrib;
   call symputx('tpmMin',min);
   call symputx('tpmQ1',q1);
   call symputx('tpmQ3',q3);
   call symputx('tpmMax',max);
run;

data export_data;
  retain transcript_id tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 tpm_bin mean_tpm;
  set tpm_data;
  if mean_tpm le &tpmQ1. then tpm_bin=1;
  else if mean_tpm le &tpmQ3. then tpm_bin=2;
  else tpm_bin=3;
  diff_tpm=(tpm_nsc1-tpm_nsc2);
  mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
  diff_log_tpm=(log_tpm_nsc1-log_tpm_nsc2);
run;

proc sort data=export_data;
   by tpm_bin;
proc corr data=export_data pearson;
  by tpm_bin;
  where mean_tpm > 0;
  var log_tpm_nsc1 log_tpm_nsc2;
run;

proc export data=export_data outfile="!MCLAB/event_analysis/analysis_output/rsem_fpkm_&xslist..csv"
   dbms=csv replace;
run;

%mend;


%concord(pb_plus_50p_apn0);
%concord(pb_plus_75p_apn0);
%concord(pb_plus_50p_apn5);
%concord(pb_plus_75p_apn5);

/*			LOW		MED		HIGH
LIST			r2	N	r2	N	r2	N
pb_plus_50p_apn0	0.19488	23765	0.76285	21306	0.93494	10664
pb_plus_50p_apn0	0.23512	17699	0.75150	17967	0.93511	8977
pb_plus_75p_apn5	0.45314	11242	0.71684	13183	0.93494	6596
pb_plus_75p_apn5	0.56775	8317	0.70333	11240	0.93577	5621

*/


