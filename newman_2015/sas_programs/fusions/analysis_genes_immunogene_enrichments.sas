/* Enrichment of fusions for T1D genes and autoimmune genes */

libname con '/home/jrbnewman/concannon/sas_data';

data gene_results_for_enrich;
   set con.results_by_fusion_genes_w_flags;
  
   if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_on_any=1;
   else flag_on_any=0;
   
   if flag_cd4_on=1 and flag_cd8_on=1 and flag_cd19_on=1 then flag_gene_on_all=1;
   else flag_gene_on_all=0;

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then do;
              if flag_gene_on_all=1 then flag_fdr05_any=0;
              else flag_fdr05_any=0;
              end;
   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fdr05_any=1;
   else flag_fdr05_any=0;


   if flag_immuno_gene=. then flag_immuno_gene=0;
   if flag_diabetes_gene=. then flag_diabetes_gene=0;
run;

ods listing;
ods html close;


/* DE enrichment - Autoimmune genes */
proc freq data=gene_results_for_enrich;
   where flag_gene_on_all=1 and flag_pseudogene=0;
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
             0 |    986 |     42 |   1028
               |   6.28 |   0.27 |   6.55
               |  95.91 |   4.09 |
               |   6.94 |   2.82 |
      ---------+--------+--------+
             1 |  13220 |   1446 |  14666
               |  84.24 |   9.21 |  93.45
               |  90.14 |   9.86 |
               |  93.06 |  97.18 |
      ---------+--------+--------+
      Total       14206     1488    15694
                  90.52     9.48   100.00 

  Statistics for Table of flag_fdr05_any by flag_immuno_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     37.3171    <.0001
    Likelihood Ratio Chi-Square    1     45.8787    <.0001
    Continuity Adj. Chi-Square     1     36.6473    <.0001
    Mantel-Haenszel Chi-Square     1     37.3147    <.0001
    Phi Coefficient                       0.0488
    Contingency Coefficient               0.0487
    Cramer's V                            0.0488


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)       986
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 15694



*/
        
/* DE enrichment - T1D genes */
proc freq data=gene_results_for_enrich;
   where flag_gene_on_all=1 and flag_pseudogene=0;
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
            0 |   1027 |      1 |   1028
              |   6.54 |   0.01 |   6.55
              |  99.90 |   0.10 |
              |   6.56 |   3.03 |
     ---------+--------+--------+
            1 |  14634 |     32 |  14666
              |  93.25 |   0.20 |  93.45
              |  99.78 |   0.22 |
              |  93.44 |  96.97 |
     ---------+--------+--------+
     Total       15661       33    15694
                 99.79     0.21   100.00

                The SAS System      16:06 Thurs

Statistics for Table of flag_fdr05_any by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      0.6694    0.4133
   Likelihood Ratio Chi-Square    1      0.8261    0.3634
   Continuity Adj. Chi-Square     1      0.2171    0.6412
   Mantel-Haenszel Chi-Square     1      0.6693    0.4133
   Phi Coefficient                       0.0065
   Contingency Coefficient               0.0065
   Cramer's V                            0.0065

    WARNING: 25% of the cells have expected counts less
             than 5. Chi-Square may not be a valid test.


 Statistics for Table of flag_fdr05_any by flag_diabetes_gene

                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      1027
              Left-sided Pr <= F          0.8933
              Right-sided Pr >= F         0.3540

              Table Probability (P)       0.2473
              Two-sided Pr <= P           0.7224

                     Sample Size = 15694

*/

/* Detection enrichment - Autoimmune genes */
proc freq data=gene_results_for_enrich;
    where flag_pseudogene=0;
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
         0 |   4607 |    174 |   4781
           |  21.05 |   0.79 |  21.84
           |  96.36 |   3.64 |
           |  22.86 |  10.03 |
  ---------+--------+--------+
         1 |  15547 |   1560 |  17107
           |  71.03 |   7.13 |  78.16
           |  90.88 |   9.12 |
           |  77.14 |  89.97 |
  ---------+--------+--------+
  Total       20154     1734    21888
              92.08     7.92   100.00

             The SAS System      16:06 Th


  Statistics for Table of flag_on_any by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    153.8143    <.0001
   Likelihood Ratio Chi-Square    1    180.3656    <.0001
   Continuity Adj. Chi-Square     1    153.0641    <.0001
   Mantel-Haenszel Chi-Square     1    153.8073    <.0001
   Phi Coefficient                       0.0838
   Contingency Coefficient               0.0835
   Cramer's V                            0.0838


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      4607
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 21888

*/
        
/* Detection enrichment - T1D genes */
proc freq data=gene_results_for_enrich;
    where flag_pseudogene=0;
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
            0 |   4780 |      1 |   4781
              |  21.84 |   0.00 |  21.84
              |  99.98 |   0.02 |
              |  21.87 |   2.94 |
     ---------+--------+--------+
            1 |  17074 |     33 |  17107
              |  78.01 |   0.15 |  78.16
              |  99.81 |   0.19 |
              |  78.13 |  97.06 |
     ---------+--------+--------+
     Total       21854       34    21888
                 99.84     0.16   100.00

                The SAS System      16:06 Thursda

 Statistics for Table of flag_on_any by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      7.1266    0.0076
   Likelihood Ratio Chi-Square    1     10.2964    0.0013
   Continuity Adj. Chi-Square     1      6.0608    0.0138
   Mantel-Haenszel Chi-Square     1      7.1263    0.0076
   Phi Coefficient                       0.0180
   Contingency Coefficient               0.0180
   Cramer's V                            0.0180


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      4780
             Left-sided Pr <= F          0.9998
             Right-sided Pr >= F         0.0024

             Table Probability (P)       0.0022
             Two-sided Pr <= P           0.0055

                    Sample Size = 21888


*/
