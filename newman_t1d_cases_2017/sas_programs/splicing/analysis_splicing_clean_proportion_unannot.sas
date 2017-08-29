/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';



/* Proportion of Unannotated events of all detected events */


data splicing_clean_w_de_flags2;
   set splicing_clean_w_de_flags;
   if flag_junction_annotated=0 and flag_intron_retention=0 then flag_junction_unannotated=1;
   else flag_junction_unannotated=0;
run;

* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1;
   tables flag_junction_unannotated*flag_immuno_gene / chisq ;
run;


/*

 Table of flag_junction_unannotated by flag_immuno_gene

          flag_junction_unannotated
                    flag_immuno_gene

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       0|       1|  Total
          ---------+--------+--------+
                 0 | 150981 |  20189 | 171170
                   |  71.93 |   9.62 |  81.55
                   |  88.21 |  11.79 |
                   |  81.50 |  81.92 |
          ---------+--------+--------+
                 1 |  34262 |   4455 |  38717
                   |  16.32 |   2.12 |  18.45
                   |  88.49 |  11.51 |
                   |  18.50 |  18.08 |
          ---------+--------+--------+
          Total      185243    24644   209887
                      88.26    11.74   100.00


 Statistics for Table of flag_junction_unannotated by flag_immuno_gene

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1      2.5296    0.1117
         Likelihood Ratio Chi-Square    1      2.5410    0.1109
         Continuity Adj. Chi-Square     1      2.5019    0.1137
         Mantel-Haenszel Chi-Square     1      2.5296    0.1117
         Phi Coefficient                      -0.0035
         Contingency Coefficient               0.0035
         Cramer's V                           -0.0035


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)    150981
                   Left-sided Pr <= F          0.0566
                   Right-sided Pr >= F         0.9454

                   Table Probability (P)       0.0020
                   Two-sided Pr <= P           0.1136

                          Sample Size = 209887

*/

* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_junction_unannotated*flag_diabetes_gene / chisq ;
run;


/*
  Table of flag_junction_unannotated by flag_diabetes_gene

            flag_junction_unannotated
                      flag_diabetes_gene

            Frequency|
            Percent  |
            Row Pct  |
            Col Pct  |       0|       1|  Total
            ---------+--------+--------+
                   0 |  19811 |    378 |  20189
                     |  80.39 |   1.53 |  81.92
                     |  98.13 |   1.87 |
                     |  81.93 |  81.47 |
            ---------+--------+--------+
                   1 |   4369 |     86 |   4455
                     |  17.73 |   0.35 |  18.08
                     |  98.07 |   1.93 |
                     |  18.07 |  18.53 |
            ---------+--------+--------+
            Total       24180      464    24644
                        98.12     1.88   100.00

                       The SAS System         11:48 Monday
Statistics for Table of flag_junction_unannotated by flag_diabetes_gene

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1      0.0667    0.7962
         Likelihood Ratio Chi-Square    1      0.0663    0.7968
         Continuity Adj. Chi-Square     1      0.0390    0.8435
         Mantel-Haenszel Chi-Square     1      0.0667    0.7962
         Phi Coefficient                       0.0016
         Contingency Coefficient               0.0016
         Cramer's V                            0.0016


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)     19811
                   Left-sided Pr <= F          0.6294
                   Right-sided Pr >= F         0.4171

                   Table Probability (P)       0.0465
                   Two-sided Pr <= P           0.8075

                          Sample Size = 24644


*/


/* IR+Unannot vs Annotated junctions */


* AUTOIMMUNE vs ALL OTHER ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1;
   tables flag_junction_annotated*flag_immuno_gene / chisq ;
run;

/*

 Table of flag_junction_annotated by flag_immuno_gene

         flag_junction_annotated
                   flag_immuno_gene

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |  70460 |   9149 |  79609
                  |  33.57 |   4.36 |  37.93
                  |  88.51 |  11.49 |
                  |  38.04 |  37.12 |
         ---------+--------+--------+
                1 | 114783 |  15495 | 130278
                  |  54.69 |   7.38 |  62.07
                  |  88.11 |  11.89 |
                  |  61.96 |  62.88 |
         ---------+--------+--------+
         Total      185243    24644   209887
                     88.26    11.74   100.00

Statistics for Table of flag_junction_annotated by flag_immuno_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      7.6819    0.0056
       Likelihood Ratio Chi-Square    1      7.7007    0.0055
       Continuity Adj. Chi-Square     1      7.6433    0.0057
       Mantel-Haenszel Chi-Square     1      7.6819    0.0056
       Phi Coefficient                       0.0060
       Contingency Coefficient               0.0060
       Cramer's V                            0.0060


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)     70460
                 Left-sided Pr <= F          0.9973
                 Right-sided Pr >= F         0.0028

                 Table Probability (P)       0.0001
                 Two-sided Pr <= P           0.0057

                        Sample Size = 209887

*/

* T1D vs AUTOIMMUNE ;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_junction_annotated*flag_diabetes_gene / chisq ;
run;


/*
  Table of flag_junction_annotated by flag_diabetes_gene

           flag_junction_annotated
                     flag_diabetes_gene

           Frequency|
           Percent  |
           Row Pct  |
           Col Pct  |       0|       1|  Total
           ---------+--------+--------+
                  0 |   8951 |    198 |   9149
                    |  36.32 |   0.80 |  37.12
                    |  97.84 |   2.16 |
                    |  37.02 |  42.67 |
           ---------+--------+--------+
                  1 |  15229 |    266 |  15495
                    |  61.80 |   1.08 |  62.88
                    |  98.28 |   1.72 |
                    |  62.98 |  57.33 |
           ---------+--------+--------+
           Total       24180      464    24644
                       98.12     1.88   100.00

                      The SAS System         11:48 Monday,

 Statistics for Table of flag_junction_annotated by flag_diabetes_gene

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1      6.2354    0.0125
         Likelihood Ratio Chi-Square    1      6.1269    0.0133
         Continuity Adj. Chi-Square     1      5.9955    0.0143
         Mantel-Haenszel Chi-Square     1      6.2352    0.0125
         Phi Coefficient                      -0.0159
         Contingency Coefficient               0.0159
         Cramer's V                           -0.0159


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)      8951
                   Left-sided Pr <= F          0.0075
                   Right-sided Pr >= F         0.9942

                   Table Probability (P)       0.0018
                   Two-sided Pr <= P           0.0133

                          Sample Size = 24644


*/



