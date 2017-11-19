ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";


/* Export RSEM FPKM data for making plots -- I also want to bin transcripts into low/low-med/med-high/high exp bins:
   bin 1 = lowest is <1
   bin 2 = lowest is between 1 and 2
   bin 3 = lowest is between 2 and 4
   bin 4 = lowest is >4
   */

%macro export(datain);

data tpm_data;
   set event.rsem_&datain.;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
   min_log_tpm=min(log_tpm_nsc1,log_tpm_nsc2);
   mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
   min_tpm=min(tpm_nsc1,tpm_nsc2);
   keep transcript_id mean_tpm tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 min_tpm min_log_tpm mean_log_tpm ;
run;


data export_data;
  set tpm_data;

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


proc export data=export_data outfile="!MCLAB/event_analysis/analysis_output/rsem_tpm_&datain..csv"
dbms=csv replace;
run;

%mend;

%export(pacbio_all);
%export(refseq_all);
%export(events_exp_any);
%export(events_exp_100perc);
%export(events_exp_75perc);
%export(events_exp_75perc_apn5);
%export(events_exp_100perc_apn5);
%export(events_exp_75perc_apn10);
%export(events_exp_100perc_apn10);
%export(events_exp_50perc);
%export(events_exp_50perc_apn5);
%export(events_exp_50perc_apn10);
%export(pacbio_apn0);
%export(pacbio_apn5);
%export(pacbio_apn10);

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

*/



