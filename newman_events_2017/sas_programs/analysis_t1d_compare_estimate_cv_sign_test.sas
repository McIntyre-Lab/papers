ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For the transcripts in common between the all and filtered list, calculate the difference in variance
   and perform the wilcoxon sign test */

data var_all;
  set eventloc.hg19_cv_all_xs;
  keep transcript_id cv_cd19 cv_cd4 cv_cd8;
  rename cv_cd4=cv_cd4_all cv_cd8=cv_cd8_all cv_cd19=cv_cd19_all;
run;

data var_filtered;
   set eventloc.hg19_cv_filtered_xs;
  keep transcript_id cv_cd19 cv_cd4 cv_cd8;
  rename cv_cd4=cv_cd4_filter cv_cd8=cv_cd8_filter cv_cd19=cv_cd19_filter;
run;

proc sort data=var_all;
   by transcript_id;
proc sort data=var_filtered;
   by transcript_id;
run;

data xs_var_all_data;
   merge  var_all (in=in1) var_filtered (in=in2);
   by transcript_id;
   if in1 and in2 ;
run;

/* Calc diff between filtered and unfiltered for each cell type, and flag -1 if lower and +1 if higher */

data calc_flag_var_diff;
   set xs_var_all_data;
   if cv_cd4_all=. or cv_cd4_filter=. then diff_cd4=.;
   else diff_cd4=cv_cd4_all-cv_cd4_filter;
   if cv_cd8_all=. or cv_cd8_filter=. then diff_cd8=.;
   else diff_cd8=cv_cd8_all-cv_cd8_filter;
   if cv_cd19_all=. or cv_cd19_filter=. then diff_cd19=.;
   else diff_cd19=cv_cd19_all-cv_cd19_filter;

   if diff_cd4=. then flag_diff_cd4=.;
      else if diff_cd4 < 0 then flag_diff_cd4=-1;
      else if diff_cd4 > 0 then flag_diff_cd4=1;
      else flag_diff_cd4=0;

   if diff_cd8=. then flag_diff_cd8=.;
      else if diff_cd8 < 0 then flag_diff_cd8=-1;
      else if diff_cd8 > 0 then flag_diff_cd8=1;
      else flag_diff_cd8=0;

   if diff_cd19=. then flag_diff_cd19=.;
      else if diff_cd19 < 0 then flag_diff_cd19=-1;
      else if diff_cd19 > 0 then flag_diff_cd19=1;
      else flag_diff_cd19=0;
run;

/* Check flags -- this will give us an idea as how big the difference will be */

proc freq data=calc_flag_var_diff;
  tables flag_diff_cd4 flag_diff_cd8 flag_diff_cd19 ;
run;

/*

1=CV is higher in complete set, -1=CV is higher in restricted set, 0=no difference 

                                            Cumulative    Cumulative
  flag_diff_cd4    Frequency     Percent     Frequency      Percent
  ------------------------------------------------------------------
             -1        3696       14.17          3696        14.17
              0          67        0.26          3763        14.42
              1       22329       85.58         26092       100.00

                        Frequency Missing = 345


                                            Cumulative    Cumulative
  flag_diff_cd8    Frequency     Percent     Frequency      Percent
  ------------------------------------------------------------------
             -1        4143       15.88          4143        15.88
              0          44        0.17          4187        16.04
              1       21910       83.96         26097       100.00

                        Frequency Missing = 340


                                              Cumulative    Cumulative
   flag_diff_cd19    Frequency     Percent     Frequency      Percent
   -------------------------------------------------------------------
               -1        3260       12.47          3260        12.47
                0          49        0.19          3309        12.66
                1       22827       87.34         26136       100.00

                         Frequency Missing = 301


*/

proc univariate data=calc_flag_var_diff;
   var diff_cd4 diff_cd8 diff_cd19 ;
run;

/* OUTPUT:
diff_CD4:
            Tests for Location: Mu0=0
 Test           -Statistic-    -----p Value------

 Student's t    t  37.47335    Pr > |t|    <.0001
 Sign           M    9316.5    Pr >= |M|   <.0001
 Signed Rank    S   1.204E8    Pr >= |S|   <.0001

diff_CD8:
            Tests for Location: Mu0=0
Test           -Statistic-    -----p Value------
 Student's t    t   38.4216    Pr > |t|    <.0001
 Sign           M    8883.5    Pr >= |M|   <.0001
 Signed Rank    S  1.1646E8    Pr >= |S|   <.0001

diff_CD19:
            Tests for Location: Mu0=0
 Test           -Statistic-    -----p Value------
Student's t    t  39.96746    Pr > |t|    <.0001
Sign           M    9783.5    Pr >= |M|   <.0001
Signed Rank    S  1.2591E8    Pr >= |S|   <.0001

*/

