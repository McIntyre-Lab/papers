/* Enrichment of genes for T1D genes and autoimmune genes */

/* Redo'ing without exlcuding pseudogenes */

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
   where flag_gene_on_all=1;
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
           0 |   6721 |     42 |   6763
             |  20.32 |   0.13 |  20.45
             |  99.38 |   0.62 |
             |  21.28 |   2.82 |
    ---------+--------+--------+
           1 |  24861 |   1446 |  26307
             |  75.18 |   4.37 |  79.55
             |  94.50 |   5.50 |
             |  78.72 |  97.18 |
    ---------+--------+--------+
    Total       31582     1488    33070
                95.50     4.50   100.00


 Statistics for Table of flag_fdr05_any by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    297.6187    <.0001
   Likelihood Ratio Chi-Square    1    425.7713    <.0001
   Continuity Adj. Chi-Square     1    296.4851    <.0001
   Mantel-Haenszel Chi-Square     1    297.6097    <.0001
   Phi Coefficient                       0.0949
   Contingency Coefficient               0.0944
   Cramer's V                            0.0949


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      6721
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 33070



*/
        
/* DE enrichment - T1D genes */
proc freq data=gene_results_for_enrich;
   where flag_gene_on_all=1;
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
             0 |   6762 |      1 |   6763
               |  20.45 |   0.00 |  20.45
               |  99.99 |   0.01 |
               |  20.47 |   3.03 |
      ---------+--------+--------+
             1 |  26275 |     32 |  26307
               |  79.45 |   0.10 |  79.55
               |  99.88 |   0.12 |
               |  79.53 |  96.97 |
      ---------+--------+--------+
      Total       33037       33    33070
                  99.90     0.10   100.00


Statistics for Table of flag_fdr05_any by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      6.1619    0.0131
   Likelihood Ratio Chi-Square    1      8.8607    0.0029
   Continuity Adj. Chi-Square     1      5.1366    0.0234
   Mantel-Haenszel Chi-Square     1      6.1617    0.0131
   Phi Coefficient                       0.0137
   Contingency Coefficient               0.0136
   Cramer's V                            0.0137


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      6762
             Left-sided Pr <= F          0.9995
             Right-sided Pr >= F         0.0050

             Table Probability (P)       0.0044
             Two-sided Pr <= P           0.0084

                    Sample Size = 33070







*/

/* Detection enrichment - Autoimmune genes */
proc freq data=gene_results_for_enrich;
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
         0 |  25181 |    174 |  25355
           |  39.20 |   0.27 |  39.47
           |  99.31 |   0.69 |
           |  40.29 |  10.03 |
  ---------+--------+--------+
         1 |  37321 |   1560 |  38881
           |  58.10 |   2.43 |  60.53
           |  95.99 |   4.01 |
           |  59.71 |  89.97 |
  ---------+--------+--------+
  Total       62502     1734    64236
              97.30     2.70   100.00

             The SAS System      16:06 Thu

Statistics for Table of flag_on_any by flag_immuno_gene

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1    646.3648    <.0001
 Likelihood Ratio Chi-Square    1    777.3029    <.0001
 Continuity Adj. Chi-Square     1    645.0991    <.0001
 Mantel-Haenszel Chi-Square     1    646.3547    <.0001
 Phi Coefficient                       0.1003
 Contingency Coefficient               0.0998
 Cramer's V                            0.1003


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     25181
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                  Sample Size = 64236





*/
        
/* Detection enrichment - T1D genes */
proc freq data=gene_results_for_enrich;
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
           0 |  25354 |      1 |  25355
             |  39.47 |   0.00 |  39.47
             | 100.00 |   0.00 |
             |  39.49 |   2.94 |
    ---------+--------+--------+
           1 |  38848 |     33 |  38881
             |  60.48 |   0.05 |  60.53
             |  99.92 |   0.08 |
             |  60.51 |  97.06 |
    ---------+--------+--------+
    Total       64202       34    64236
                99.95     0.05   100.00

               The SAS System      16:06 Thursday, Dec

             The FREQ Procedure


Statistics for Table of flag_on_any by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     19.0009    <.0001
  Likelihood Ratio Chi-Square    1     25.9821    <.0001
  Continuity Adj. Chi-Square     1     17.5019    <.0001
  Mantel-Haenszel Chi-Square     1     19.0006    <.0001
  Phi Coefficient                       0.0172
  Contingency Coefficient               0.0172
  Cramer's V                            0.0172


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     25354
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 64236



*/
