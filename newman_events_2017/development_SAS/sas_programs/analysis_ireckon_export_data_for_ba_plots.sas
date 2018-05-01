ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Export iReckon data for making BA plots. I am going to calculate the data needed for plots here then export
   to python for plotting */

data NSC1;
   set event.ireckon_results_nsc;
   where sample_id="NSC1";
   length target_id $200.;
   target_id=catx("|",gene_id,transcript_id);
   keep target_id rpkm;
   rename target_id=transcript_id rpkm=rpkm_nsc1;
run;

data NSC2;
   set event.ireckon_results_nsc;
   where sample_id="NSC2";
   length target_id $200.;
   target_id=catx("|",gene_id,transcript_id);
   keep target_id rpkm;
   rename target_id=transcript_id rpkm=rpkm_nsc2;
run;

proc sort data=NSC1 nodup;
  by transcript_id;
proc sort data=NSC2 nodup;
  by transcript_id;
run;


proc means data=NSC1 noprint;
  by transcript_id;
  var rpkm_nsc1 ;
  output out=NSC1_2(drop=_TYPE_ _FREQ_) max=;
run;

proc means data=NSC2 noprint;
  by  transcript_id;
  var rpkm_nsc2 ;
  output out=NSC2_2(drop=_TYPE_ _FREQ_) max=;
run;




data calc_mean_diff;
  merge NSC1_2 (in=in1) NSC2_2 (in=in2);
  by transcript_id;
  if in1 and in2;
  *if rpkm_nsc1=. then rpkm_nsc1=0;
  *if rpkm_nsc2=. then rpkm_nsc2=0;
  log_rpkm_nsc1=log(rpkm_nsc1+1);
  log_rpkm_nsc2=log(rpkm_nsc2+1);
  mean_rpkm=(rpkm_nsc1+rpkm_nsc2)/2;
  min_log_rpkm=min(log_rpkm_nsc1,log_rpkm_nsc2);
  mean_log_rpkm=(log_rpkm_nsc1+log_rpkm_nsc2)/2;
  min_rpkm=min(rpkm_nsc1,rpkm_nsc2);
  keep transcript_id mean_rpkm rpkm_nsc1 rpkm_nsc2 log_rpkm_nsc1 log_rpkm_nsc2 min_rpkm min_log_rpkm mean_log_rpkm ;
run;

/* look at logRPKM distribution for non-zero samples */
proc univariate data=calc_mean_diff noprint;
   where mean_rpkm > 0;
   var min_log_rpkm;
   output out=rpkm_distrib min=min q1=q1 q3=q3 p90=p90 max=max;
run;

proc print data=rpkm_distrib;
run;




/*

   max        p90         q3         q1              min

 8.37088    1.61964    0.65689    0.070142    .000010067


BIN	DEFINITION	logRPKM RANGE	Approximate Range
0	meanTPM=0	0				0
1	0->Q1		0-0.07			0 - 0.1
2	Q1->Q3		0.07-0.65	0.1 - 0.6
3	Q3->90%		0.65-1.61	0.5 - 1.6
4	90%->100%	1.6-8.37	1.6 - Max

*/


data export_data;
  set calc_mean_diff;

  /* bin 4 levels */
  if mean_log_rpkm=0 then rpkm_bin=0;
  else if min_log_rpkm < 0.1 then rpkm_bin=1;
  else if min_log_rpkm < 0.6 then rpkm_bin=2;
  else if min_log_rpkm < 1.6 then rpkm_bin=3;
  else rpkm_bin=4;

  diff_rpkm=abs(rpkm_nsc1-rpkm_nsc2);
  diff_log_rpkm=abs(log_rpkm_nsc1-log_rpkm_nsc2);
  diff_over_mean_log_rpkm=diff_log_rpkm/mean_log_rpkm;
run;


proc sort data=export_data;
   by rpkm_bin;
proc corr data=export_data out=pearson_log pearson;
  by rpkm_bin;
  var log_rpkm_nsc1 log_rpkm_nsc2;
run;   

proc sort data=export_data;
   by rpkm_bin;
proc corr data=export_data out=pearson_raw pearson;
  by rpkm_bin;
  var rpkm_nsc1 rpkm_nsc2;
run;   

data r_log;
  length dataset $30.;
  length stat $30.;
  set pearson_log;
  where _TYPE_="N" or _NAME_="log_rpkm_nsc1";
  dataset="iReckon";
  stat=catt("bin_",rpkm_bin,"_",_TYPE_);
  keep dataset stat log_rpkm_nsc2;
  rename log_rpkm_nsc2=log_value;
run;

data r_raw;
  length dataset $30.;
  length stat $30.;
  set pearson_raw;
  where _TYPE_="N" or _NAME_="rpkm_nsc1";
  dataset="iReckon";
  stat=catt("bin_",rpkm_bin,"_",_TYPE_);
  keep dataset stat rpkm_nsc2;
  rename rpkm_nsc2=raw_value;
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
proc transpose data=corr_data out=r2_sbys_ireckon;
   by dataset;
   id stat;
   var r2_raw r2_log;
run;


proc sort data=export_data;
   by descending rpkm_bin;
run;


proc export data=export_data outfile="!MCLAB/event_analysis/analysis_output/ireckon_rpkm_complete_common.csv"
dbms=csv replace;
run;

proc print data=r2_sbys_ireckon;
     where _NAME_ = "r2_log";
run;

proc print data=r2_sbys_ireckon;
     where _NAME_ = "r2_raw";
run;


/*

          bin_1_             bin_2_             bin_3_            bin_4_
 _NAME_    CORR    bin_1_N    CORR    bin_2_N    CORR   bin_3_N    CORR   bin_4_N

 r2_log  0.013875    801    0.013193    996    0.10419    414    0.82158    251

                                                                  bin_4_
 _NAME_ bin_1_CORR bin_1_N bin_2_CORR bin_2_N bin_3_CORR bin_3_N   CORR  bin_4_N

 r2_raw .000391206   801   .000232970   996   .001225139   414   0.90628   251

*/

