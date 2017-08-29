/* Need to determine what percentage of exons are differentially expressed (DE)

Do for: all genes, autoimmune genes, T1D genes
Then we want the proportion of DE autoimmune among all genes, and proportion of DE T1D among autoimmune 

Going to use three "DE" definitions:

1. DE exons that are expressed among all tissues
2. DE exons expressed among all tissues PLUS tissue-specific exons
3. DE exons expressed among all tissues PLUS tissue-specific exons PLUS two-tissue exons

We can then decide which to report. Option 3 is probably the most complete definition of DE because:
	- it includes exons that are DE between all three tissues
	- it includes exons that are specific to a single tissue (by definition, these will be DE)
	- it includes exons that are expressed in two tissues only (by definition, these will also be DE)

*/

/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Prep data */

data fus_results_for_enrich;
   set con.results_by_fusion_w_flags;

  /* Scenario 1 - only DE fusion expressed in all three cell types */

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then flag_fusion_de=.;
   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fusion_de=1;
   else flag_fusion_de=0;

  /* Scenario 2 - only DE fusion expressed in all three cell types PLUS fusions expressed in only one cell type */

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then do;
        if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_fusion_dd_1=1;
        else if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_fusion_dd_1=1;
        else if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_fusion_dd_1=1;
        else flag_fusion_dd_1=.;
        end;

   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fusion_dd_1=1;
   else flag_fusion_dd_1=0;

  /* Scenario 3 - only DE fusion expressed in all three cell types PLUS fusions expressed in only one or two cell types */

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then do;
     if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_fusion_dd_2=1;
     else flag_fusion_dd_2=.;
     end;
   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fusion_dd_2=1;
   else flag_fusion_dd_2=0;
   

   if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_any_on=1;
   else flag_any_on=0;
   if flag_immuno_gene=. then flag_immuno_gene=0;
   if flag_diabetes_gene=. then flag_diabetes_gene=0;
run;

ods listing;
ods html close;

/* Get DE and DD counts */

*all genes - DE, exp all 3 tissues;

proc freq data=fus_results_for_enrich;
   where flag_any_on=1;
   tables flag_fusion_de;
run;

/* 

                                            Cumulative    Cumulative
 flag_fusion_de    Frequency     Percent     Frequency      Percent
 -------------------------------------------------------------------
              0       38554       23.55         38554        23.55
              1      125159       76.45        163713       100.00

                      Frequency Missing = 27676

*/


*all genes - DE + tissue-specific ;

proc freq data=fus_results_for_enrich;
   where flag_any_on=1;
   tables flag_fusion_dd_1;
run;

/* 

                                             Cumulative    Cumulative
flag_fusion_dd_1    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       38554       23.17         38554        23.17
               1      127868       76.83        166422       100.00

                      Frequency Missing = 24967

*/


*all genes - DE + two-tissue fusions ;

proc freq data=fus_results_for_enrich;
   where flag_any_on=1;
   tables flag_fusion_dd_2;
run;

/* 

                                                Cumulative    Cumulative
   flag_fusion_dd_2    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0       38554       20.14         38554        20.14
                  1      152835       79.86        191389       100.00

*/



*Autoimmune genes - DE, exp all 3 tissues;

proc freq data=fus_results_for_enrich;
   where flag_immuno_gene=1 and flag_any_on=1;
   tables flag_fusion_de;
run;

/* 

                                             Cumulative    Cumulative
  flag_fusion_de    Frequency     Percent     Frequency      Percent
  -------------------------------------------------------------------
               0        2771       18.91          2771        18.91
               1       11885       81.09         14656       100.00

                        Frequency Missing = 1906

*/


*Autoimmune genes - DE + tissue-specific ;

proc freq data=fus_results_for_enrich;
   where flag_immuno_gene=1 and flag_any_on=1;
   tables flag_fusion_dd_1;
run;

/* 

                                               Cumulative    Cumulative
  flag_fusion_dd_1    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        2771       18.69          2771        18.69
                 1       12052       81.31         14823       100.00

                         Frequency Missing = 1739

*/


*Autoimmune genes - DE + two-tissue fusions ;

proc freq data=fus_results_for_enrich;
   where flag_immuno_gene=1 and flag_any_on=1;
   tables flag_fusion_dd_2;
run;

/* 


                                               Cumulative    Cumulative
  flag_fusion_dd_2    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        2771       16.73          2771        16.73
                 1       13791       83.27         16562       100.00
*/



*T1D genes - DE, exp all 3 tissues;

proc freq data=fus_results_for_enrich;
   where flag_diabetes_gene=1 and flag_any_on=1;
   tables flag_fusion_de;
run;

/* 


                                            Cumulative    Cumulative
 flag_fusion_de    Frequency     Percent     Frequency      Percent
 -------------------------------------------------------------------
              0          52       19.85            52        19.85
              1         210       80.15           262       100.00

                        Frequency Missing = 50

*/


*T1D genes - DE + tissue-specific ;

proc freq data=fus_results_for_enrich;
   where flag_diabetes_gene=1 and flag_any_on=1;
   tables flag_fusion_dd_1;
run;

/* 

                                               Cumulative    Cumulative
  flag_fusion_dd_1    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0          52       19.77            52        19.77
                 1         211       80.23           263       100.00

                          Frequency Missing = 49

*/


*T1D genes - DE + two-tissue fusions ;

proc freq data=fus_results_for_enrich;
   where flag_diabetes_gene=1 and flag_any_on=1;
   tables flag_fusion_dd_2;
run;

/* 

                                              Cumulative    Cumulative
 flag_fusion_dd_2    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0          52       16.67            52        16.67
                1         260       83.33           312       100.00

*/



/* Summary so far:

If only DE exons expressed in all three tissues: 76.45% all exons, 81.09% autoimmune exons, 80.15% T1D exons
If only DE exons expressed in all three tissues PLUS tissue-specific exons: 76.83% all exons, 81.31% autoimmune exons, 80.23% T1D exons
If only DE exons expressed in all three tissues PLUS two-tissue exons: 79.86% all exons, 83.27% autoimmune exons, 83.33% T1D exons

*/


/* What is the proportion of autoimmune exons among DE/DD exons? */

* only DE exons expressed in all three tissues ;

proc freq data=fus_results_for_enrich;
   where flag_fusion_all_on0=1;
   tables flag_fusion_de*flag_immuno_gene / chisq;
run;

/*
Table of flag_fusion_de by flag_immuno_gene

    flag_fusion_de     flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  35783 |   2771 |  38554
             |  21.86 |   1.69 |  23.55
             |  92.81 |   7.19 |
             |  24.01 |  18.91 |
    ---------+--------+--------+
           1 | 113274 |  11885 | 125159
             |  69.19 |   7.26 |  76.45
             |  90.50 |   9.50 |
             |  75.99 |  81.09 |
    ---------+--------+--------+
    Total      149057    14656   163713
                91.05     8.95   100.00

Statistics for Table of flag_fusion_de by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1    192.7277    <.0001
  Likelihood Ratio Chi-Square    1    201.6428    <.0001
  Continuity Adj. Chi-Square     1    192.4445    <.0001
  Mantel-Haenszel Chi-Square     1    192.7265    <.0001
  Phi Coefficient                       0.0343
  Contingency Coefficient               0.0343
  Cramer's V                            0.0343


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     35783
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 163713


*/

* only DE exons expressed in all three tissues PLUS tissue-specific exons;

proc freq data=fus_results_for_enrich;
   where flag_any_on=1;
   tables flag_fusion_dd_1*flag_immuno_gene / chisq;
run;

/*
 Table of flag_fusion_dd_1 by flag_immuno_ge

      flag_fusion_dd_1
                flag_immuno_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  35783 |   2771 |  38554
               |  21.50 |   1.67 |  23.17
               |  92.81 |   7.19 |
               |  23.60 |  18.69 |
      ---------+--------+--------+
             1 | 115816 |  12052 | 127868
               |  69.59 |   7.24 |  76.83
               |  90.57 |   9.43 |
               |  76.40 |  81.31 |
      ---------+--------+--------+
      Total      151599    14823   166422
                  91.09     8.91   100.00

Statistics for Table of flag_fusion_dd_1 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    182.8684    <.0001
   Likelihood Ratio Chi-Square    1    191.2074    <.0001
   Continuity Adj. Chi-Square     1    182.5927    <.0001
   Mantel-Haenszel Chi-Square     1    182.8673    <.0001
   Phi Coefficient                       0.0331
   Contingency Coefficient               0.0331
   Cramer's V                            0.0331

Statistics for Table of flag_fusion_dd_1 by flag_immuno_gene

                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     35783
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

               Effective Sample Size = 166422
                 Frequency Missing = 24967

*/

* only DE exons expressed in all three tissues PLUS two-tissue exons;

proc freq data=fus_results_for_enrich;
   where flag_any_on=1;
   tables flag_fusion_dd_2*flag_immuno_gene / chisq;
run;

/*

 Table of flag_fusion_dd_2 by flag_immuno_gene

      flag_fusion_dd_2
                flag_immuno_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  35783 |   2771 |  38554
               |  18.70 |   1.45 |  20.14
               |  92.81 |   7.19 |
               |  20.47 |  16.73 |
      ---------+--------+--------+
             1 | 139044 |  13791 | 152835
               |  72.65 |   7.21 |  79.86
               |  90.98 |   9.02 |
               |  79.53 |  83.27 |
      ---------+--------+--------+
      Total      174827    16562   191389
                  91.35     8.65   100.00


Statistics for Table of flag_fusion_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    131.3098    <.0001
   Likelihood Ratio Chi-Square    1    136.8485    <.0001
   Continuity Adj. Chi-Square     1    131.0776    <.0001
   Mantel-Haenszel Chi-Square     1    131.3091    <.0001
   Phi Coefficient                       0.0262
   Contingency Coefficient               0.0262
   Cramer's V                            0.0262


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     35783
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 191389

*/


/* What is the proportion of T1D exons among DE/DD autoimmune exons? */


* only DE exons expressed in all three tissues ;

proc freq data=fus_results_for_enrich;
   where flag_fusion_all_on0=1 and flag_immuno_gene=1;
   tables flag_fusion_de*flag_diabetes_gene / chisq;
run;

/*

 Table of flag_fusion_de by flag_diabetes_gene

      flag_fusion_de
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   2719 |     52 |   2771
               |  18.55 |   0.35 |  18.91
               |  98.12 |   1.88 |
               |  18.89 |  19.85 |
      ---------+--------+--------+
             1 |  11675 |    210 |  11885
               |  79.66 |   1.43 |  81.09
               |  98.23 |   1.77 |
               |  81.11 |  80.15 |
      ---------+--------+--------+
      Total       14394      262    14656
                  98.21     1.79   100.00


 Statistics for Table of flag_fusion_de by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      0.1539    0.6949
    Likelihood Ratio Chi-Square    1      0.1520    0.6966
    Continuity Adj. Chi-Square     1      0.0978    0.7545
    Mantel-Haenszel Chi-Square     1      0.1539    0.6949
    Phi Coefficient                      -0.0032
    Contingency Coefficient               0.0032
    Cramer's V                           -0.0032


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      2719
              Left-sided Pr <= F          0.3719
              Right-sided Pr >= F         0.6858

              Table Probability (P)       0.0577
              Two-sided Pr <= P           0.6907

                     Sample Size = 14656



*/

* only DE exons expressed in all three tissues PLUS tissue-specific exons;

proc freq data=fus_results_for_enrich;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_fusion_dd_1*flag_diabetes_gene / chisq;
run;

/*
 Table of flag_fusion_dd_1 by flag_diabetes_gene

       flag_fusion_dd_1
                 flag_diabetes_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |   2719 |     52 |   2771
                |  18.34 |   0.35 |  18.69
                |  98.12 |   1.88 |
                |  18.67 |  19.77 |
       ---------+--------+--------+
              1 |  11841 |    211 |  12052
                |  79.88 |   1.42 |  81.31
                |  98.25 |   1.75 |
                |  81.33 |  80.23 |
       ---------+--------+--------+
       Total       14560      263    14823
                   98.23     1.77   100.00

 Statistics for Table of flag_fusion_dd_1 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1      0.2047    0.6510
     Likelihood Ratio Chi-Square    1      0.2018    0.6533
     Continuity Adj. Chi-Square     1      0.1389    0.7094
     Mantel-Haenszel Chi-Square     1      0.2047    0.6510
     Phi Coefficient                      -0.0037
     Contingency Coefficient               0.0037
     Cramer's V                           -0.0037


           Fisher's Exact Test
    ----------------------------------
    Cell (1,1) Frequency (F)      2719
    Left-sided Pr <= F          0.3496
    Right-sided Pr >= F         0.7066

    Table Probability (P)       0.0563
    Two-sided Pr <= P           0.6326

      Effective Sample Size = 14823
         Frequency Missing = 1739

*/

* only DE exons expressed in all three tissues PLUS two-tissue exons;

proc freq data=fus_results_for_enrich;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_fusion_dd_2*flag_diabetes_gene / chisq;
run;

/*
  Table of flag_fusion_dd_2 by flag_diabetes_gene

        flag_fusion_dd_2
                  flag_diabetes_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |   2719 |     52 |   2771
                 |  16.42 |   0.31 |  16.73
                 |  98.12 |   1.88 |
                 |  16.73 |  16.67 |
        ---------+--------+--------+
               1 |  13531 |    260 |  13791
                 |  81.70 |   1.57 |  83.27
                 |  98.11 |   1.89 |
                 |  83.27 |  83.33 |
        ---------+--------+--------+
        Total       16250      312    16562
                    98.12     1.88   100.00

                   The SAS System        14:16 Monday, D

 Statistics for Table of flag_fusion_dd_2 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1      0.0009    0.9755
     Likelihood Ratio Chi-Square    1      0.0009    0.9754
     Continuity Adj. Chi-Square     1      0.0000    1.0000
     Mantel-Haenszel Chi-Square     1      0.0009    0.9755
     Phi Coefficient                       0.0002
     Contingency Coefficient               0.0002
     Cramer's V                            0.0002


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)      2719
               Left-sided Pr <= F          0.5363
               Right-sided Pr >= F         0.5248

               Table Probability (P)       0.0611
               Two-sided Pr <= P           1.0000

                      Sample Size = 16562

*/


/* Summary:

If only DE exons expressed in all three tissues:
76.45% all exons, 81.09% autoimmune exons, 80.15% T1D exons
81.09% of autoimmune exons vs. 75.99% of non-autoimmune exons
80.15% of T1D exons vs. 81.11% of non-T1D autoimmune exons


If only DE exons expressed in all three tissues PLUS tissue-specific exons:
76.83% all exons, 81.31% autoimmune exons, 80.23% T1D exons
81.31% of autoimmune exons vs. 76.40% of non-autoimmune exons
80.23% of T1D exons vs. 81.33% of non-T1D autoimmune exons


If only DE exons expressed in all three tissues PLUS two-tissue exons:
79.86% all exons DE, 83.27% autoimmune exons DE, 83.33% T1D exons DE
83.27% of autoimmune exons vs. 79.53% of non-autoimmune exons
83.33% of T1D exons vs. 83.27% of non-T1D autoimmune exons

*/

