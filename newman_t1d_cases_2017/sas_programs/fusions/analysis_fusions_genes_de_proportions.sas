/* Need to determine what percentage of genes are differentially expressed (DE) using Scenario 3

Do for: all genes, autoimmune genes, T1D genes
Then we want the proportion of DE autoimmune among all genes, and proportion of DE T1D among autoimmune 

DE definition: DE exons expressed among all tissues PLUS tissue-specific exons PLUS two-tissue exons

*/

/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';

data gene_results_for_enrich;
   set con.results_by_fusion_genes_w_flags;
  
   if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_any_on=1;
   else flag_any_on=0;
   
   if flag_cd4_on=1 and flag_cd8_on=1 and flag_cd19_on=1 then flag_gene_on_all=1;
   else flag_gene_on_all=0;

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then do;
              if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_gene_de=1;
               else flag_gene_de=.; *if gene is completely OFF, then it shouldn't have a DE flag;
              end;
   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_gene_de=1;
   else flag_gene_de=0;

   if flag_immuno_gene=. then flag_immuno_gene=0;
   if flag_diabetes_gene=. then flag_diabetes_gene=0;
run;

ods listing;
ods html close;


/* Get DE counts */

* all genes;
proc freq data=gene_results_for_enrich;
   where flag_any_on=1 and flag_pseudogene=0;
   tables flag_gene_de;
run;

/*
                                          Cumulative    Cumulative
 flag_gene_de    Frequency     Percent     Frequency      Percent
 -----------------------------------------------------------------
            0        1014        5.93          1014         5.93
            1       16093       94.07         17107       100.00


*/

* autoimmune genes;
proc freq data=gene_results_for_enrich;
   where flag_any_on=1 and flag_immuno_gene=1 and flag_pseudogene=0;
   tables flag_gene_de;
run;

/*
                                      Cumulative    Cumulative
g_gene_de    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------
        0          41        2.63            41         2.63
        1        1519       97.37          1560       100.00
*/

* T1D genes;
proc freq data=gene_results_for_enrich;
   where flag_any_on=1 and flag_diabetes_gene=1 and flag_pseudogene=0;
   tables flag_gene_de;
run;

/*
                                        Cumulative    Cumulative
lag_gene_de    Frequency     Percent     Frequency      Percent
----------------------------------------------------------------
          0           1        3.03             1         3.03
          1          32       96.97            33       100.00


*/

/* What proportion of DE genes are autoimmune genes? */

proc freq data=gene_results_for_enrich;
   where flag_any_on=1 and flag_pseudogene=0;
   tables flag_gene_de*flag_immuno_gene / chisq;
run;

/*
 Table of flag_gene_de by flag_immuno_gene

    flag_gene_de     flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |    973 |     41 |   1014
             |   5.69 |   0.24 |   5.93
             |  95.96 |   4.04 |
             |   6.26 |   2.63 |
    ---------+--------+--------+
           1 |  14574 |   1519 |  16093
             |  85.19 |   8.88 |  94.07
             |  90.56 |   9.44 |
             |  93.74 |  97.37 |
    ---------+--------+--------+
    Total       15547     1560    17107
                90.88     9.12   100.00

 Statistics for Table of flag_gene_de by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     33.5073    <.0001
  Likelihood Ratio Chi-Square    1     41.0335    <.0001
  Continuity Adj. Chi-Square     1     32.8595    <.0001
  Mantel-Haenszel Chi-Square     1     33.5054    <.0001
  Phi Coefficient                       0.0443
  Contingency Coefficient               0.0442
  Cramer's V                            0.0443


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)       973
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 17107
*/

/* What proportion of DE autoimmune genes are T1D genes? */

proc freq data=gene_results_for_enrich;
   where flag_any_on=1 and flag_immuno_gene=1 and flag_pseudogene=0;
   tables flag_gene_de*flag_diabetes_gene / chisq;
run;

/*
 
   Table of flag_gene_de by flag_diabetes_gene

       flag_gene_de     flag_diabetes_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |     40 |      1 |     41
                |   2.56 |   0.06 |   2.63
                |  97.56 |   2.44 |
                |   2.62 |   3.03 |
       ---------+--------+--------+
              1 |   1487 |     32 |   1519
                |  95.32 |   2.05 |  97.37
                |  97.89 |   2.11 |
                |  97.38 |  96.97 |
       ---------+--------+--------+
       Total        1527       33     1560
                   97.88     2.12   100.00

 Statistics for Table of flag_gene_de by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      0.0213    0.8840
   Likelihood Ratio Chi-Square    1      0.0203    0.8866
   Continuity Adj. Chi-Square     1      0.0000    1.0000
   Mantel-Haenszel Chi-Square     1      0.0213    0.8840
   Phi Coefficient                      -0.0037
   Contingency Coefficient               0.0037
   Cramer's V                           -0.0037

    WARNING: 25% of the cells have expected counts less
             than 5. Chi-Square may not be a valid test.

 Statistics for Table of flag_gene_de by flag_diabetes_gene

                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)        40
             Left-sided Pr <= F          0.5886
             Right-sided Pr >= F         0.7857

             Table Probability (P)       0.3743
             Two-sided Pr <= P           0.5886

                     Sample Size = 1560



*/

