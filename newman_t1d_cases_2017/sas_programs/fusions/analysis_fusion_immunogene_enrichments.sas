/* Enrichment of fusions for T1D genes and autoimmune genes */

libname con '/home/jrbnewman/concannon/sas_data';

data fus_results_for_enrich;
   set con.results_by_fusion_w_flags;
   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then flag_fdr05_any=.;
   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fdr05_any=1;
   else flag_fdr05_any=0;
   
   if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_on_any=1;
   else flag_on_any=0;
   if flag_immuno_gene=. then flag_immuno_gene=0;
   if flag_diabetes_gene=. then flag_diabetes_gene=0;
run;

ods listing;
ods html close;


/* DE enrichment - Autoimmune genes */
proc freq data=fus_results_for_enrich;
   where flag_fusion_all_on0=1;
   tables flag_fdr05_any*flag_immuno_gene / chisq;
run;

/*
 Table of flag_fdr05_any by flag_immuno_gene

     flag_fdr05_any     flag_immuno_gene

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

                The SAS System    07:45 Thur

 Statistics for Table of flag_fdr05_any by flag_immuno_gene

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
        
/* DE enrichment - T1D genes */
proc freq data=fus_results_for_enrich;
   where flag_fusion_all_on0=1;
   tables flag_fdr05_any*flag_diabetes_gene / chisq;
run;


/*
 Table of flag_fdr05_any by flag_diabetes_gene

      flag_fdr05_any
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  38502 |     52 |  38554
               |  23.52 |   0.03 |  23.55
               |  99.87 |   0.13 |
               |  23.56 |  19.85 |
      ---------+--------+--------+
             1 | 124949 |    210 | 125159
               |  76.32 |   0.13 |  76.45
               |  99.83 |   0.17 |
               |  76.44 |  80.15 |
      ---------+--------+--------+
      Total      163451      262   163713
                  99.84     0.16   100.00

Statistics for Table of flag_fdr05_any by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      1.9980    0.1575
   Likelihood Ratio Chi-Square    1      2.0775    0.1495
   Continuity Adj. Chi-Square     1      1.7974    0.1800
   Mantel-Haenszel Chi-Square     1      1.9980    0.1575
   Phi Coefficient                       0.0035
   Contingency Coefficient               0.0035
   Cramer's V                            0.0035


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     38502
             Left-sided Pr <= F          0.9337
             Right-sided Pr >= F         0.0882

             Table Probability (P)       0.0218
             Two-sided Pr <= P           0.1666

                    Sample Size = 163713



*/

/* Detection enrichment - Autoimmune genes */
proc freq data=fus_results_for_enrich;
   tables flag_on_any*flag_immuno_gene / chisq;
run;


/*
  Table of flag_on_any by flag_immuno_gene

    flag_on_any     flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 | 133163 |   5005 | 138168
             |  40.41 |   1.52 |  41.93
             |  96.38 |   3.62 |
             |  43.24 |  23.21 |
    ---------+--------+--------+
           1 | 174827 |  16562 | 191389
             |  53.05 |   5.03 |  58.07
             |  91.35 |   8.65 |
             |  56.76 |  76.79 |
    ---------+--------+--------+
    Total      307990    21567   329557
                93.46     6.54   100.00


 Statistics for Table of flag_on_any by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1   3320.9935    <.0001
  Likelihood Ratio Chi-Square    1   3550.7816    <.0001
  Continuity Adj. Chi-Square     1   3320.1709    <.0001
  Mantel-Haenszel Chi-Square     1   3320.9834    <.0001
  Phi Coefficient                       0.1004
  Contingency Coefficient               0.0999
  Cramer's V                            0.1004


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    133163
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 329557


*/
        
/* Detection enrichment - T1D genes */
proc freq data=fus_results_for_enrich;
   tables flag_on_any*flag_diabetes_gene / chisq;
run;



/*
 Table of flag_on_any by flag_diabetes_gene

    flag_on_any     flag_diabetes_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 | 138106 |     62 | 138168
             |  41.91 |   0.02 |  41.93
             |  99.96 |   0.04 |
             |  41.95 |  16.58 |
    ---------+--------+--------+
           1 | 191077 |    312 | 191389
             |  57.98 |   0.09 |  58.07
             |  99.84 |   0.16 |
             |  58.05 |  83.42 |
    ---------+--------+--------+
    Total      329183      374   329557
                99.89     0.11   100.00

               The SAS System    07:45 Thursday

  Statistics for Table of flag_on_any by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     98.8059    <.0001
    Likelihood Ratio Chi-Square    1    111.0654    <.0001
    Continuity Adj. Chi-Square     1     97.7664    <.0001
    Mantel-Haenszel Chi-Square     1     98.8056    <.0001
    Phi Coefficient                       0.0173
    Contingency Coefficient               0.0173
    Cramer's V                            0.0173


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)    138106
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 329557

*/
