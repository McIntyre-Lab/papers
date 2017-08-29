/*** AUTOIMMUNE AND T1D SPLICING EVENTS - PROPORTION OF EVENTS ***/

/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';



/* Proportion of exon skipping of all detected events */

* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_exonskip*flag_immuno_gene / chisq;
run;

/*
   Table of flag_exonskip by flag_immuno_gene

     flag_exonskip     flag_immuno_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 | 144124 |  19682 | 163806
              |  68.67 |   9.38 |  78.04
              |  87.98 |  12.02 |
              |  77.80 |  79.87 |
     ---------+--------+--------+
            1 |  41119 |   4962 |  46081
              |  19.59 |   2.36 |  21.96
              |  89.23 |  10.77 |
              |  22.20 |  20.13 |
     ---------+--------+--------+
     Total      185243    24644   209887
                 88.26    11.74   100.00


 Statistics for Table of flag_exonskip by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     54.0035    <.0001
   Likelihood Ratio Chi-Square    1     54.9794    <.0001
   Continuity Adj. Chi-Square     1     53.8832    <.0001
   Mantel-Haenszel Chi-Square     1     54.0032    <.0001
   Phi Coefficient                      -0.0160
   Contingency Coefficient               0.0160
   Cramer's V                           -0.0160


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    144124
             Left-sided Pr <= F          <.0001
             Right-sided Pr >= F         1.0000

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 209887


*/

* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_exonskip*flag_diabetes_gene / chisq;
run;


/*
   Table of flag_exonskip by flag_diabetes_gene

      flag_exonskip
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  19329 |    353 |  19682
               |  78.43 |   1.43 |  79.87
               |  98.21 |   1.79 |
               |  79.94 |  76.08 |
      ---------+--------+--------+
             1 |   4851 |    111 |   4962
               |  19.68 |   0.45 |  20.13
               |  97.76 |   2.24 |
               |  20.06 |  23.92 |
      ---------+--------+--------+
      Total       24180      464    24644
                  98.12     1.88   100.00

                 The SAS System         11:48 Mo

Statistics for Table of flag_exonskip by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      4.2191    0.0400
   Likelihood Ratio Chi-Square    1      4.0429    0.0444
   Continuity Adj. Chi-Square     1      3.9824    0.0460
   Mantel-Haenszel Chi-Square     1      4.2189    0.0400
   Phi Coefficient                       0.0131
   Contingency Coefficient               0.0131
   Cramer's V                            0.0131


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     19329
             Left-sided Pr <= F          0.9811
             Right-sided Pr >= F         0.0247

             Table Probability (P)       0.0058
             Two-sided Pr <= P           0.0467

                    Sample Size = 24644

*/


/* Proportion of IR events of all detected events */

* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_intron_retention*flag_immuno_gene / chisq;
run;


/*
  Table of flag_intron_retention by flag_immuno_gene

        flag_intron_retention
                  flag_immuno_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 | 149045 |  19950 | 168995
                 |  71.01 |   9.51 |  80.52
                 |  88.19 |  11.81 |
                 |  80.46 |  80.95 |
        ---------+--------+--------+
               1 |  36198 |   4694 |  40892
                 |  17.25 |   2.24 |  19.48
                 |  88.52 |  11.48 |
                 |  19.54 |  19.05 |
        ---------+--------+--------+
        Total      185243    24644   209887
                    88.26    11.74   100.00

                   The SAS System         11:48 Mond

                 The FREQ Procedure

 Statistics for Table of flag_intron_retention by flag_immuno_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      3.3780    0.0661
       Likelihood Ratio Chi-Square    1      3.3947    0.0654
       Continuity Adj. Chi-Square     1      3.3466    0.0673
       Mantel-Haenszel Chi-Square     1      3.3779    0.0661
       Phi Coefficient                      -0.0040
       Contingency Coefficient               0.0040
       Cramer's V                           -0.0040


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)    149045
                 Left-sided Pr <= F          0.0334
                 Right-sided Pr >= F         0.9678

                 Table Probability (P)       0.0013
                 Two-sided Pr <= P           0.0670

                        Sample Size = 209887


*/


* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_intron_retention*flag_diabetes_gene / chisq;
run;


/*
 
Table of flag_intron_retention by flag_diabetes_gene

        flag_intron_retention
                  flag_diabetes_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |  19598 |    352 |  19950
                 |  79.52 |   1.43 |  80.95
                 |  98.24 |   1.76 |
                 |  81.05 |  75.86 |
        ---------+--------+--------+
               1 |   4582 |    112 |   4694
                 |  18.59 |   0.45 |  19.05
                 |  97.61 |   2.39 |
                 |  18.95 |  24.14 |
        ---------+--------+--------+
        Total       24180      464    24644
                    98.12     1.88   100.00

                   The SAS System         11:48 Monday,

Statistics for Table of flag_intron_retention by flag_diabetes_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      7.9481    0.0048
       Likelihood Ratio Chi-Square    1      7.4829    0.0062
       Continuity Adj. Chi-Square     1      7.6152    0.0058
       Mantel-Haenszel Chi-Square     1      7.9478    0.0048
       Phi Coefficient                       0.0180
       Contingency Coefficient               0.0180
       Cramer's V                            0.0180


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)     19598
                 Left-sided Pr <= F          0.9975
                 Right-sided Pr >= F         0.0036

                 Table Probability (P)       0.0010
                 Two-sided Pr <= P           0.0060

                        Sample Size = 24644



*/


