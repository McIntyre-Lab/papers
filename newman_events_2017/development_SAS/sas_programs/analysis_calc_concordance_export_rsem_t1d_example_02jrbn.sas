ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Data for example concordance and BA plots for T1D data */

%macro concord(datain,dataout);

data sample1;
   set eventloc.&datain.;
   where library="SL30324";
   keep tpm transcript_id;
   rename tpm=tpm_sample1;
run;

data sample2;
   set eventloc.&datain.;
   where library="SL30327";
   keep tpm transcript_id;
   rename tpm=tpm_sample2;
run;

proc sort data=sample1;
  by transcript_id;
proc sort data=sample2;
  by transcript_id;
run;

data merged;
  merge sample1 (in=in1) sample2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data tpm_data;
  set merged;
  log_tpm_sample1=log(tpm_sample1+1);
  log_tpm_sample2=log(tpm_sample2+1);
  mean_tpm=(tpm_sample1+tpm_sample2)/2;
  min_log_tpm=min(log_tpm_sample1,log_tpm_sample2);
  mean_log_tpm=(log_tpm_sample1+log_tpm_sample2)/2;
  min_tpm=min(tpm_sample1,tpm_sample2);
  keep transcript_id mean_tpm tpm_sample1 tpm_sample2 log_tpm_sample1 log_tpm_sample2 min_tpm min_log_tpm mean_log_tpm ;
run;

data export_data;
  set tpm_data;

  /* bin 4 levels */
  if mean_log_tpm=0 then tpm_bin=0;
  else if min_log_tpm < 0.5 then tpm_bin=1;
  else if min_log_tpm < 2 then tpm_bin=2;
  else if min_log_tpm < 4 then tpm_bin=3;
  else tpm_bin=4;

  diff_tpm=abs(tpm_sample1-tpm_sample2);
  diff_log_tpm=abs(log_tpm_sample1-log_tpm_sample2);
run;

/* CV calculation */
proc sort data=export_data;
  by tpm_bin;
proc means data=export_data noprint;
  by tpm_bin;
  var mean_log_tpm;
  output out=var_stats_by_bin
         var(mean_log_tpm)=var_mean_log  cv(mean_log_tpm)=cv_mean_log;
run;

proc print data=var_stats_by_bin;
run;

/* Agreement test */


data agreement_data;
  set export_data;
  if tpm_sample1=0 then flag_sample1_tpm_gt0=0;
  else flag_sample1_tpm_gt0=1;
  if tpm_sample2=0 then flag_sample2_tpm_gt0=0;
  else flag_sample2_tpm_gt0=1;
run;

proc freq data=agreement_data ;
   where flag_sample1_tpm_gt0 ne 0 and flag_sample2_tpm_gt0 ne 0;
   tables flag_sample1_tpm_gt0*flag_sample2_tpm_gt0 / out=counts;
   test agree;
   output out=kappa_stats AGREE;
run;


proc sort data=export_data;
   by descending tpm_bin;
run;

data export_data2;
  retain transcript_id;
  set export_data;
  diff_over_mean_log_tpm=0;
run;

proc export data=export_data2 outfile="!MCLAB/event_analysis/analysis_output/rsem_t1d_concordance_example_&dataout..csv"
   dbms=csv replace;
run;

%mend;

%concord(hg19_rsem_all_xscripts,all);
%concord(hg19_rsem_75perc_apn5_xscripts,reduced);


/* CV:
FULL:
                                        var_      cv_mean_
Obs    tpm_bin    _TYPE_    _FREQ_    mean_log       log

 1        0          0      304244     0.00000       .
 2        1          0      231058     0.03088     91.8347
 3        2          0       47033     0.18383     36.9941
 4        3          0       12107     0.28273     18.8017
 5        4          0        1645     1.02576     19.9563

 flag_sample1_tpm_gt0
           flag_sample2_tpm_gt0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 304244 |  33059 | 337303
          |  51.04 |   5.55 |  56.59
          |  90.20 |   9.80 |
          |  78.26 |  15.95 |
 ---------+--------+--------+
        1 |  84514 | 174270 | 258784
          |  14.18 |  29.24 |  43.41
          |  32.66 |  67.34 |
          |  21.74 |  84.05 |
 ---------+--------+--------+
 Total      388758   207329   596087
             65.22    34.78   100.00

            McNemar's Test
      ---------------------------
      Statistic (S)    22518.9204
      DF                        1
      Pr > S               <.0001


        Simple Kappa Coefficient
    --------------------------------
    Kappa                     0.5890
    ASE                       0.0011
    95% Lower Conf Limit      0.5870
    95% Upper Conf Limit      0.5911

         Test of H0: Kappa = 0

    ASE under H0              0.0013
    Z                       462.3177
    One-sided Pr >  Z         <.0001
    Two-sided Pr > |Z|        <.0001

          Sample Size = 596087


REDUCED:

                                         var_      cv_mean_
 Obs    tpm_bin    _TYPE_    _FREQ_    mean_log       log

  1        0          0       3729      0.00000       .
  2        1          0       6343      0.06579     65.9713
  3        2          0       8641      0.20049     31.8352
  4        3          0       6115      0.29933     18.5377
  5        4          0       1609      1.26554     21.5457



   flag_sample1_tpm_gt0
             flag_sample2_tpm_gt0

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |   3729 |   1379 |   5108
            |  14.11 |   5.22 |  19.32
            |  73.00 |  27.00 |
            |  63.92 |   6.69 |
   ---------+--------+--------+
          1 |   2105 |  19224 |  21329
            |   7.96 |  72.72 |  80.68
            |   9.87 |  90.13 |
            |  36.08 |  93.31 |
   ---------+--------+--------+
   Total        5834    20603    26437
               22.07    77.93   100.00


              McNemar's Test
         -------------------------
         Statistic (S)    151.2847
         DF                      1
         Pr > S             <.0001


          Simple Kappa Coefficient
      --------------------------------
      Kappa                     0.5990
      ASE                       0.0061
      95% Lower Conf Limit      0.5870
      95% Upper Conf Limit      0.6109

           Test of H0: Kappa = 0

      ASE under H0              0.0061
      Z                        97.7307
      One-sided Pr >  Z         <.0001
      Two-sided Pr > |Z|        <.0001

            Sample Size = 26437



*/

