ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For the transcripts in common between the all and filtered list, calculate the difference in variance
   and perform the wilcoxon sign test */


data xs2keep;
   set event.hg19_flag_xscripts_w_unique;
   if perc_features_dtct_cd4 ge 0.75
   or perc_features_dtct_cd8 ge 0.75
   or perc_features_dtct_cd19 ge 0.75;
   keep transcript_id;
run;


data var_all;
  set eventloc.hg19_variance_all_xs;
  keep transcript_id var_cd19 var_cd4 var_cd8;
  rename var_cd4=var_cd4_all var_cd8=var_cd8_all var_cd19=var_cd19_all;
run;

data var_filtered;
   set eventloc.hg19_variance_filtered_xs;
  keep transcript_id var_cd19 var_cd4 var_cd8;
  rename var_cd4=var_cd4_filter var_cd8=var_cd8_filter var_cd19=var_cd19_filter;
run;

proc sort data=xs2keep;
   by transcript_id;
proc sort data=var_all;
   by transcript_id;
proc sort data=var_filtered;
   by transcript_id;
run;

data xs_var_all_data;
   merge xs2keep (in=in1) var_all (in=in2) var_filtered (in=in3);
   by transcript_id;
   if in1 and in2 and in3;
run;

/* Calc diff between filtered and unfiltered for each cell type, and flag -1 if lower and +1 if higher */

data calc_flag_var_diff;
   set xs_var_all_data;
   diff_cd4=var_cd4_all-var_cd4_filter;
   diff_cd8=var_cd8_all-var_cd8_filter;
   diff_cd19=var_cd19_all-var_cd19_filter;
   if diff_cd4 < 0 then flag_diff_cd4=-1;
      else if diff_cd4 > 0 then flag_diff_cd4=1;
      else flag_diff_cd4=0;
   if diff_cd8 < 0 then flag_diff_cd8=-1;
      else if diff_cd8 > 0 then flag_diff_cd8=1;
      else flag_diff_cd8=0;
   if diff_cd19 < 0 then flag_diff_cd19=-1;
      else if diff_cd19 > 0 then flag_diff_cd19=1;
      else flag_diff_cd19=0;
run;

/* Check flags -- this will give us an idea as how big the difference will be */

proc freq data=calc_flag_var_diff;
  tables flag_diff_cd4 flag_diff_cd8 flag_diff_cd19 ;
run;

/*

                                             Cumulative    Cumulative
   flag_diff_cd4    Frequency     Percent     Frequency      Percent
   ------------------------------------------------------------------
              -1       23900       91.27         23900        91.27
               0         663        2.53         24563        93.80
               1        1623        6.20         26186       100.00


                                             Cumulative    Cumulative
   flag_diff_cd8    Frequency     Percent     Frequency      Percent
   ------------------------------------------------------------------
              -1       23351       89.17         23351        89.17
               0         651        2.49         24002        91.66
               1        2184        8.34         26186       100.00


                                             Cumulative    Cumulative
  flag_diff_cd19    Frequency     Percent     Frequency      Percent
  -------------------------------------------------------------------
              -1       23763       90.75         23763        90.75
               0         685        2.62         24448        93.36
               1        1738        6.64         26186       100.00

*/

proc univariate data=calc_flag_var_diff;
   var diff_cd4 diff_cd8 diff_cd19 ;
run;

/* OUTPUT:
diff_CD4:
 Test           -Statistic-    -----p Value------

 Sign           M  -11138.5    Pr >= |M|   <.0001
 Signed Rank    S   -1.51E8    Pr >= |S|   <.0001


diff_CD8:
Test           -Statistic-    -----p Value------

Sign           M  -10583.5    Pr >= |M|   <.0001
Signed Rank    S  -1.446E8    Pr >= |S|   <.0001

diff_CD19:
    Test           -Statistic-    -----p Value------

    Sign           M  -11012.5    Pr >= |M|   <.0001
    Signed Rank    S  -1.485E8    Pr >= |S|   <.0001

*/

