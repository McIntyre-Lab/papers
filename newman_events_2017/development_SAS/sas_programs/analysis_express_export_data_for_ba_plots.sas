ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Export eXpress data for making BA plots. I am going to calculate the data needed for plots here then export
   to python for plotting

   Put in macro so I can iterate through the data */

%macro export(datain);

data NSC1;
   set event.xprs_NSC1_&datain.;
   keep target_id tpm;
   rename target_id=transcript_id tpm=tpm_nsc1;
run;

data NSC2;
   set event.xprs_NSC2_&datain.;
   keep target_id tpm;
   rename target_id=transcript_id tpm=tpm_nsc2;
run;

proc sort data=NSC1;
  by transcript_id;
proc sort data=NSC2;
  by transcript_id;
run;

data calc_mean_diff;
  merge NSC1 (in=in1) NSC2 (in=in2);
  by transcript_id;
  if in1 and in2;
  log_tpm_nsc1=log(tpm_nsc1+1);
  log_tpm_nsc2=log(tpm_nsc2+1);
  mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
  min_log_tpm=min(log_tpm_nsc1,log_tpm_nsc2);
  mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
  min_tpm=min(tpm_nsc1,tpm_nsc2);
  keep transcript_id mean_tpm tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 min_tpm min_log_tpm mean_log_tpm ;
run;


data export_data;
  set calc_mean_diff;

  /* bin 4 levels */
  if mean_log_tpm=0 then tpm_bin=0;
  else if min_log_tpm < 0.5 then tpm_bin=1;
  else if min_log_tpm < 2 then tpm_bin=2;
  else if min_log_tpm < 4 then tpm_bin=3;
  else tpm_bin=4;

  diff_tpm=abs(tpm_nsc1-tpm_nsc2);
  diff_log_tpm=abs(log_tpm_nsc1-log_tpm_nsc2);
  diff_over_mean_log_tpm=diff_log_tpm/mean_log_tpm;
run;


proc sort data=export_data;
   by tpm_bin;
proc corr data=export_data out=pearson_log pearson;
  by tpm_bin;
  var log_tpm_nsc1 log_tpm_nsc2;
run;   

proc sort data=export_data;
   by tpm_bin;
proc corr data=export_data out=pearson_raw pearson;
  by tpm_bin;
  var tpm_nsc1 tpm_nsc2;
run;   


data r_log;
  length dataset $30.;
  length stat $30.;
  set pearson_log;
  where _TYPE_="N" or _NAME_="log_tpm_nsc1";
  dataset="&datain.";
  stat=catt("bin_",tpm_bin,"_",_TYPE_);
  keep dataset stat log_tpm_nsc2;
  rename log_tpm_nsc2=log_value;
run;

data r_raw;
  length dataset $30.;
  length stat $30.;
  set pearson_raw;
  where _TYPE_="N" or _NAME_="tpm_nsc1";
  dataset="&datain.";
  stat=catt("bin_",tpm_bin,"_",_TYPE_);
  keep dataset stat tpm_nsc2;
  rename tpm_nsc2=raw_value;
run;

proc sort data=r_log;
  by dataset stat;
proc sort data=r_raw;
  by dataset stat;
run;

data corr_data;
  merge r_log r_raw;
  by dataset stat;
  if index(stat,"_CORR") > 0 then do;
     r2_raw=raw_value ** 2;
     r2_log=log_value ** 2;
     end;
  else do; 
    r2_raw=raw_value;
    r2_log=log_value;
  end;

run;


proc sort data=corr_data;
   by dataset stat;
proc transpose data=corr_data out=r2_sbys_&datain.;
   by dataset;
   id stat;
   var r2_raw r2_log;
run;


proc sort data=export_data;
   by descending tpm_bin;
run;


proc export data=export_data outfile="!MCLAB/event_analysis/analysis_output/express_tpm_&datain..csv"
dbms=csv replace;
run;

%mend;

%export(refseq_all);
%export(events_100prc_apn0);
%export(events_75prc_apn5);


data all_r2_data;
  set r2_: ;
run;

proc print data=all_r2_data;
     where _NAME_ = "r2_log";
run;

proc print data=all_r2_data;
     where _NAME_ = "r2_raw";
run;

/*
dataset		   _NAME_  bin1_r  bin1_n   bin2_r bin2_n  bin3_r bin3_n  bin4_r bin4_n bin0_corr bin0_n
events_100prc_apn0 r2_log 0.11708     146  0.68249   5553 0.78161  5974  0.95080   2067         .      .
events_75prc_apn5  r2_log 0.13414    1735  0.62611   5128 0.75683  5648  0.95037   2106         .    117
refseq_all         r2_log 0.61976   64319  0.75345  22314 0.80734  8349  0.94640   2360         .  31289

dataset		   _NAME_   bin1_r  bin1_n   bin2_r bin2_n  bin3_r bin3_n  bin4_r bin4_n bin0_corr bin0_n
events_100prc_apn0 r2_raw 0.062679     146  0.37165   5553 0.07321   5974 0.95513   2067         .      .
events_75prc_apn5  r2_raw 0.012072    1735  0.02889   5128 0.08523   5648 0.95073   2106         .    117
refseq_all         r2_raw 0.019534   64319  0.11577  22314 0.17920   8349 0.95005   2360         .  31289
*/

