/*** AUTOIMMUNE AND T1D SPLICING EVENTS ENRICHMENTS - EXONSKIP and INTRON RETENTION ***/


/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';



/* Proportion of detected ES events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1;
   tables flag_any_on*flag_immuno_gene / chisq;
run;


/*
Table of flag_any_on by flag_immuno_gene

  flag_any_on     flag_immuno_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |4514982 | 566184 |5081166
           |  88.06 |  11.04 |  99.10
           |  88.86 |  11.14 |
           |  99.10 |  99.13 |
  ---------+--------+--------+
         1 |  41119 |   4962 |  46081
           |   0.80 |   0.10 |   0.90
           |  89.23 |  10.77 |
           |   0.90 |   0.87 |
  ---------+--------+--------+
  Total     4556101   571146  5127247
              88.86    11.14   100.00


 Statistics for Table of flag_any_on by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      6.4808    0.0109
  Likelihood Ratio Chi-Square    1      6.5443    0.0105
  Continuity Adj. Chi-Square     1      6.4430    0.0111
  Mantel-Haenszel Chi-Square     1      6.4808    0.0109
  Phi Coefficient                      -0.0011
  Contingency Coefficient               0.0011
  Cramer's V                           -0.0011


                   Fisher's Exact Test
           -----------------------------------
           Cell (1,1) Frequency (F)    4514982
           Left-sided Pr <= F           0.0054
           Right-sided Pr >= F          0.9948

           Table Probability (P)        0.0002
           Two-sided Pr <= P            0.0107

                  Sample Size = 5127247





*/

/* Proportion of DE ES events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
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
             0 |  14857 |   1553 |  16410
               |  32.24 |   3.37 |  35.61
               |  90.54 |   9.46 |
               |  36.13 |  31.30 |
      ---------+--------+--------+
             1 |  26262 |   3409 |  29671
               |  56.99 |   7.40 |  64.39
               |  88.51 |  11.49 |
               |  63.87 |  68.70 |
      ---------+--------+--------+
      Total       41119     4962    46081
                  89.23    10.77   100.00


Statistics for Table of flag_event_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     45.1198    <.0001
   Likelihood Ratio Chi-Square    1     45.9146    <.0001
   Continuity Adj. Chi-Square     1     44.9092    <.0001
   Mantel-Haenszel Chi-Square     1     45.1188    <.0001
   Phi Coefficient                       0.0313
   Contingency Coefficient               0.0313
   Cramer's V                            0.0313


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     14857
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 46081

*/





/* Proportion of detected IR events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1;
   tables flag_any_on*flag_immuno_gene / chisq;
run;


/*

Table of flag_any_on by flag_immuno_gene

  flag_any_on     flag_immuno_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 168108 |  16180 | 184288
           |  74.65 |   7.19 |  81.84
           |  91.22 |   8.78 |
           |  82.28 |  77.51 |
  ---------+--------+--------+
         1 |  36198 |   4694 |  40892
           |  16.08 |   2.08 |  18.16
           |  88.52 |  11.48 |
           |  17.72 |  22.49 |
  ---------+--------+--------+
  Total      204306    20874   225180
              90.73     9.27   100.00

Statistics for Table of flag_any_on by flag_immuno_gene

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1    289.9178    <.0001
 Likelihood Ratio Chi-Square    1    275.6400    <.0001
 Continuity Adj. Chi-Square     1    289.5970    <.0001
 Mantel-Haenszel Chi-Square     1    289.9165    <.0001
 Phi Coefficient                       0.0359
 Contingency Coefficient               0.0359
 Cramer's V                            0.0359


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)    168108
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                  Sample Size = 225180


*/

/* Proportion of DE IR events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
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
             0 |  15350 |   1906 |  17256
               |  37.54 |   4.66 |  42.20
               |  88.95 |  11.05 |
               |  42.41 |  40.61 |
      ---------+--------+--------+
             1 |  20848 |   2788 |  23636
               |  50.98 |   6.82 |  57.80
               |  88.20 |  11.80 |
               |  57.59 |  59.39 |
      ---------+--------+--------+
      Total       36198     4694    40892
                  88.52    11.48   100.00

Statistics for Table of flag_event_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      5.5233    0.0188
   Likelihood Ratio Chi-Square    1      5.5406    0.0186
   Continuity Adj. Chi-Square     1      5.4498    0.0196
   Mantel-Haenszel Chi-Square     1      5.5232    0.0188
   Phi Coefficient                       0.0116
   Contingency Coefficient               0.0116
   Cramer's V                            0.0116


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     15350
             Left-sided Pr <= F          0.9911
             Right-sided Pr >= F         0.0097

             Table Probability (P)       0.0008
             Two-sided Pr <= P           0.0193

                    Sample Size = 40892

*/





/* Proportion of detected autoimmune ES events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_immuno_gene=1;
   tables flag_any_on*flag_diabetes_gene / chisq;
run;


/*
     Table of flag_any_on by flag_diabetes_gene

        flag_any_on     flag_diabetes_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 | 559380 |   6804 | 566184
                 |  97.94 |   1.19 |  99.13
                 |  98.80 |   1.20 |
                 |  99.14 |  98.39 |
        ---------+--------+--------+
               1 |   4851 |    111 |   4962
                 |   0.85 |   0.02 |   0.87
                 |  97.76 |   2.24 |
                 |   0.86 |   1.61 |
        ---------+--------+--------+
        Total      564231     6915   571146
                    98.79     1.21   100.00

                   The SAS System         11:48 Mo


Statistics for Table of flag_any_on by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     44.0779    <.0001
  Likelihood Ratio Chi-Square    1     35.3568    <.0001
  Continuity Adj. Chi-Square     1     43.2166    <.0001
  Mantel-Haenszel Chi-Square     1     44.0778    <.0001
  Phi Coefficient                       0.0088
  Contingency Coefficient               0.0088
  Cramer's V                            0.0088


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    559380
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 571146

*/

/* Proportion of autoimmune DE ES events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_diabetes_gene / chisq;
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
             0 |   1531 |     22 |   1553
               |  30.85 |   0.44 |  31.30
               |  98.58 |   1.42 |
               |  31.56 |  19.82 |
      ---------+--------+--------+
             1 |   3320 |     89 |   3409
               |  66.91 |   1.79 |  68.70
               |  97.39 |   2.61 |
               |  68.44 |  80.18 |
      ---------+--------+--------+
      Total        4851      111     4962
                  97.76     2.24   100.00


  Statistics for Table of flag_event_dd_2 by flag_diabetes_gene

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1      6.9566    0.0084
      Likelihood Ratio Chi-Square    1      7.5535    0.0060
      Continuity Adj. Chi-Square     1      6.4213    0.0113
      Mantel-Haenszel Chi-Square     1      6.9552    0.0084
      Phi Coefficient                       0.0374
      Contingency Coefficient               0.0374
      Cramer's V                            0.0374


                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)      1531
                Left-sided Pr <= F          0.9978
                Right-sided Pr >= F         0.0044

                Table Probability (P)       0.0022
                Two-sided Pr <= P           0.0093

                        Sample Size = 4962


*/





/* Proportion of detected autoimmune IR events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_immuno_gene=1;
   tables flag_any_on*flag_diabetes_gene / chisq;
run;


/*
Table of flag_any_on by flag_diabetes_gene

   flag_any_on     flag_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  15919 |    261 |  16180
            |  76.26 |   1.25 |  77.51
            |  98.39 |   1.61 |
            |  77.65 |  69.97 |
   ---------+--------+--------+
          1 |   4582 |    112 |   4694
            |  21.95 |   0.54 |  22.49
            |  97.61 |   2.39 |
            |  22.35 |  30.03 |
   ---------+--------+--------+
   Total       20501      373    20874
               98.21     1.79   100.00


Statistics for Table of flag_any_on by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     12.3855    0.0004
  Likelihood Ratio Chi-Square    1     11.5729    0.0007
  Continuity Adj. Chi-Square     1     11.9490    0.0005
  Mantel-Haenszel Chi-Square     1     12.3849    0.0004
  Phi Coefficient                       0.0244
  Contingency Coefficient               0.0244
  Cramer's V                            0.0244


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     15919
            Left-sided Pr <= F          0.9997
            Right-sided Pr >= F         0.0004

            Table Probability (P)       0.0001
            Two-sided Pr <= P           0.0007

                   Sample Size = 20874


*/

/* Proportion of autoimmune DE IR events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_diabetes_gene / chisq;
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
             0 |   1883 |     23 |   1906
               |  40.12 |   0.49 |  40.61
               |  98.79 |   1.21 |
               |  41.10 |  20.54 |
      ---------+--------+--------+
             1 |   2699 |     89 |   2788
               |  57.50 |   1.90 |  59.39
               |  96.81 |   3.19 |
               |  58.90 |  79.46 |
      ---------+--------+--------+
      Total        4582      112     4694
                  97.61     2.39   100.00

 Statistics for Table of flag_event_dd_2 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     19.1620    <.0001
     Likelihood Ratio Chi-Square    1     20.9127    <.0001
     Continuity Adj. Chi-Square     1     18.3190    <.0001
     Mantel-Haenszel Chi-Square     1     19.1579    <.0001
     Phi Coefficient                       0.0639
     Contingency Coefficient               0.0638
     Cramer's V                            0.0639


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)      1883
               Left-sided Pr <= F          1.0000
               Right-sided Pr >= F         <.0001

               Table Probability (P)       <.0001
               Two-sided Pr <= P           <.0001

                       Sample Size = 4694

 

*/


