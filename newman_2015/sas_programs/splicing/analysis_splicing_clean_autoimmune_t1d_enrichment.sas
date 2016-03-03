/*** AUTOIMMUNE AND T1D SPLICING EVENTS ENRICHMENTS - DETECTED AND DE ***/

/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';


/* Now need to look at proportions:

autoimmune among all detected
diabetes among autoimmune detected

autoimmune among de/dd detected
diabetes among de/dd autoimmune detected

*/

/* Autoimmune among all detected */

proc freq data=splicing_clean_w_de_flags;
   tables flag_any_on*flag_immuno_gene / chisq ;
run;

/*
 Table of flag_any_on by flag_immuno_gene

   flag_any_on     flag_immuno_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |5025413 | 621226 |5646639
            |  85.81 |  10.61 |  96.42
            |  89.00 |  11.00 |
            |  96.44 |  96.18 |
   ---------+--------+--------+
          1 | 185243 |  24644 | 209887
            |   3.16 |   0.42 |   3.58
            |  88.26 |  11.74 |
            |   3.56 |   3.82 |
   ---------+--------+--------+
   Total     5210656   645870  5856526
               88.97    11.03   100.00

 Statistics for Table of flag_any_on by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1    112.8960    <.0001
  Likelihood Ratio Chi-Square    1    110.9080    <.0001
  Continuity Adj. Chi-Square     1    112.8206    <.0001
  Mantel-Haenszel Chi-Square     1    112.8960    <.0001
  Phi Coefficient                       0.0044
  Contingency Coefficient               0.0044
  Cramer's V                            0.0044


                   Fisher's Exact Test
           -----------------------------------
           Cell (1,1) Frequency (F)    5025413
           Left-sided Pr <= F           1.0000
           Right-sided Pr >= F          <.0001

           Table Probability (P)        <.0001
           Two-sided Pr <= P            <.0001

                  Sample Size = 5856526

*/


/* T1D among autoimmune detected */

proc freq data=splicing_clean_w_de_flags;
   where flag_immuno_gene=1;
   tables flag_any_on*flag_diabetes_gene / chisq ;
run;

/*
Table of flag_any_on by flag_diabetes_gene

   flag_any_on     flag_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 | 613513 |   7713 | 621226
            |  94.99 |   1.19 |  96.18
            |  98.76 |   1.24 |
            |  96.21 |  94.33 |
   ---------+--------+--------+
          1 |  24180 |    464 |  24644
            |   3.74 |   0.07 |   3.82
            |  98.12 |   1.88 |
            |   3.79 |   5.67 |
   ---------+--------+--------+
   Total      637693     8177   645870
               98.73     1.27   100.00


Statistics for Table of flag_any_on by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     77.9711    <.0001
  Likelihood Ratio Chi-Square    1     68.2477    <.0001
  Continuity Adj. Chi-Square     1     77.4589    <.0001
  Mantel-Haenszel Chi-Square     1     77.9709    <.0001
  Phi Coefficient                       0.0110
  Contingency Coefficient               0.0110
  Cramer's V                            0.0110


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    613513
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 645870

*/


/* Autoimmune among all DE */


proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq ;
run;

/*
   Table of flag_event_dd_2 by flag_immuno_gene

       flag_event_dd_2
                 flag_immuno_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |  68492 |   8343 |  76835
                |  32.63 |   3.97 |  36.61
                |  89.14 |  10.86 |
                |  36.97 |  33.85 |
       ---------+--------+--------+
              1 | 116751 |  16301 | 133052
                |  55.63 |   7.77 |  63.39
                |  87.75 |  12.25 |
                |  63.03 |  66.15 |
       ---------+--------+--------+
       Total      185243    24644   209887
                   88.26    11.74   100.00

                  The SAS System         11:48 Mon

Statistics for Table of flag_event_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     91.2394    <.0001
   Likelihood Ratio Chi-Square    1     92.1370    <.0001
   Continuity Adj. Chi-Square     1     91.1050    <.0001
   Mantel-Haenszel Chi-Square     1     91.2390    <.0001
   Phi Coefficient                       0.0208
   Contingency Coefficient               0.0208
   Cramer's V                            0.0208


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     68492
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 209887


*/



/* Diabetes among autoimmune DE */


proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_diabetes_gene / chisq ;
run;

/*
Table of flag_event_dd_2 by flag_diabetes_gene

     flag_event_dd_2
               flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |   8255 |     88 |   8343
              |  33.50 |   0.36 |  33.85
              |  98.95 |   1.05 |
              |  34.14 |  18.97 |
     ---------+--------+--------+
            1 |  15925 |    376 |  16301
              |  64.62 |   1.53 |  66.15
              |  97.69 |   2.31 |
              |  65.86 |  81.03 |
     ---------+--------+--------+
     Total       24180      464    24644
                 98.12     1.88   100.00

                The SAS System         11:48 Mon


Statistics for Table of flag_event_dd_2 by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     46.8128    <.0001
    Likelihood Ratio Chi-Square    1     51.5640    <.0001
    Continuity Adj. Chi-Square     1     46.1376    <.0001
    Mantel-Haenszel Chi-Square     1     46.8109    <.0001
    Phi Coefficient                       0.0436
    Contingency Coefficient               0.0435
    Cramer's V                            0.0436


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      8255
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 24644

*/


