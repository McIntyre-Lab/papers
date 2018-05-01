ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For transcripts with unique features, set unique detection threshold using &percUniqOn. and drop transcripts with less than this threshold of unique features expressed. */

/* Set threshold here. I will make this into a command line variable later */
*%let percUniqOn=0.25;
*%let percUniqOn=0.5;
*%let percUniqOn=0.75;
%let percUniqOn=1;

/* %let percUniqOn=sysget(percUniqOn) */

/* Proportion (x100) will be appended to output dataset name */

%global percSuff;
%let percSuff = %SYSEVALF(&percUniqOn * 100);
%put &percSuff.;

/* Filter transcripts based on the proportion of unique features on */

data flag_xscripts_lt_perc;
   set event.flag_xscripts_w_unique;
   if flag_xscript_has_unique=1 then do;
       if perc_unique_features_dtct ge &percUniqOn. then flag_keep_xscript=1;
       else flag_keep_xscript=0;
       end;
   else do;
       if flag_xscript_has_dtct_features=1 then flag_keep_xscript=1;
       else flag_keep_xscript=0;
   end;
run;

 
/* Count number of transcripts retained */

proc freq data=flag_xscripts_lt_perc;
   tables flag_xscript_has_unique*flag_keep_xscript;
run;

  
/* Make permenant */

data event.xscripts_w_uniq_dtct_gt_&percSuff.;
   set flag_xscripts_lt_perc;
run;


/* Export transcript lists */

data xs_for_export;
  set event.xscripts_w_uniq_dtct_gt_&percSuff.;
  keep transcript_id flag_xscript_has_dtct_features flag_xscript_has_unique
       flag_xscript_has_unique_dtct flag_keep_xscript;
  run;

proc export data=xs_for_export outfile="!MCLAB/event_analysis/analysis_output/event_analysis_transcripts_kept_&percSuff.perc_unique.csv" dbms=csv replace;
run;



/*  Counts output:

For 25%:
 flag_xscript_has_unique
           flag_keep_xscript

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   2199 |  48052 |  50251
         |   2.10 |  45.81 |  47.91
         |   4.38 |  95.62 |
         |   7.93 |  62.27 |
---------+--------+--------+
       1 |  25518 |  29114 |  54632
         |  24.33 |  27.76 |  52.09
         |  46.71 |  53.29 |
         |  92.07 |  37.73 |
---------+--------+--------+
Total       27717    77166   104883
            26.43    73.57   100.00

48052 transcripts (no unique features) retained
2199 transcripts (no unique features) dropped (i.e. nothing detected!)
29114 transcripts (with unique features) retained
25518 transcripts (with unique features) dropped

For 50%:
flag_xscript_has_unique
          flag_keep_xscript

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   2199 |  48052 |  50251
         |   2.10 |  45.81 |  47.91
         |   4.38 |  95.62 |
         |   7.21 |  64.60 |
---------+--------+--------+
       1 |  28305 |  26327 |  54632
         |  26.99 |  25.10 |  52.09
         |  51.81 |  48.19 |
         |  92.79 |  35.40 |
---------+--------+--------+
Total       30504    74379   104883
            29.08    70.92   100.00

48052 transcripts (no unique features) retained
2199 transcripts (no unique features) dropped (i.e. nothing detected!)
26327 transcripts (with unique features) retained
28305 transcripts (with unique features) dropped


For 75%:

 flag_xscript_has_unique
           flag_keep_xscript

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   2199 |  48052 |  50251
          |   2.10 |  45.81 |  47.91
          |   4.38 |  95.62 |
          |   5.68 |  72.63 |
 ---------+--------+--------+
        1 |  36525 |  18107 |  54632
          |  34.82 |  17.26 |  52.09
          |  66.86 |  33.14 |
          |  94.32 |  27.37 |
 ---------+--------+--------+
 Total       38724    66159   104883
             36.92    63.08   100.00



48052 transcripts (no unique features) retained
2199 transcripts (no unique features) dropped (i.e. nothing detected!)
18107 transcripts (with unique features) retained
36525 transcripts (with unique features) dropped


For 100%:
flag_xscript_has_unique
          flag_keep_xscript

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   2199 |  48052 |  50251
         |   2.10 |  45.81 |  47.91
         |   4.38 |  95.62 |
         |   5.54 |  73.73 |
---------+--------+--------+
       1 |  37511 |  17121 |  54632
         |  35.76 |  16.32 |  52.09
         |  68.66 |  31.34 |
         |  94.46 |  26.27 |
---------+--------+--------+
Total       39710    65173   104883
            37.86    62.14   100.00

48052 transcripts (no unique features) retained
2199 transcripts (no unique features) dropped (i.e. nothing detected!)
17121 transcripts (with unique features) retained
37511 transcripts (with unique features) dropped
*/

