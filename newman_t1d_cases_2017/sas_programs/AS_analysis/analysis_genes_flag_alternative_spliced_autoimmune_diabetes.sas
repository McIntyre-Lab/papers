/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Calculate the proportion of genes with differential splicing for:
   (1) Autoimmune genes
   (2) T1D genes */

data diff_splicing;
   set con.gene_diff_splicing_summary;
run;


data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=diff_splicing;
   by gene_id;
proc sort data=ai nodup;
   by gene_id;
proc sort data=t1d nodup;
   by gene_id;
run;

data diff_splicing_immuno;
   merge diff_splicing (in=in1) ai (in=in2) t1d (in=in3);
   by gene_id;
   if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
   if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in1 then output;
run;


proc freq data=diff_splicing_immuno;
   where flag_immuno_gene=1;
   tables flag_gene_any_dd;
run;

proc freq data=diff_splicing_immuno;
   tables flag_gene_any_dd*flag_immuno_gene / chisq;
run;


/* AUTOIMMUNE GENES
                                              Cumulative    Cumulative
 flag_gene_any_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0         381       22.94           381        22.94
                1        1280       77.06          1661       100.00

                         Frequency Missing = 29

 Table of flag_gene_any_dd by flag_immuno_gene

      flag_gene_any_dd
                flag_immuno_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |  28580 |    381 |  28961
               |  64.98 |   0.87 |  65.84
               |  98.68 |   1.32 |
               |  67.53 |  22.94 |
      ---------+--------+--------+
             1 |  13744 |   1280 |  15024
               |  31.25 |   2.91 |  34.16
               |  91.48 |   8.52 |
               |  32.47 |  77.06 |
      ---------+--------+--------+
      Total       42324     1661    43985
                  96.22     3.78   100.00

            Frequency Missing = 3078

 Statistics for Table of flag_gene_any_dd by flag_immuno_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1   1412.9006    <.0001
    Likelihood Ratio Chi-Square    1   1333.2005    <.0001
    Continuity Adj. Chi-Square     1   1410.9187    <.0001
    Mantel-Haenszel Chi-Square     1   1412.8685    <.0001
    Phi Coefficient                       0.1792
    Contingency Coefficient               0.1764
    Cramer's V                            0.1792


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     28580
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                Effective Sample Size = 43985
                   Frequency Missing = 3078

*/



proc freq data=diff_splicing_immuno;
   where flag_immunobase_diabetes_gene=1;
   tables flag_gene_any_dd;
run;

proc freq data=diff_splicing_immuno;
   where flag_immuno_gene=1;
   tables flag_gene_any_dd*flag_immunobase_diabetes_gene / chisq;
run;


/* T1D GENES

                                              Cumulative    Cumulative
 flag_gene_any_dd    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0         105       26.38           105        26.38
                1         293       73.62           398       100.00

                          Frequency Missing = 7

T1D vs AUTOIMMUNE GENES

  flag_gene_any_dd
            flag_immunobase_diabetes_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |    276 |    105 |    381
           |  16.62 |   6.32 |  22.94
           |  72.44 |  27.56 |
           |  21.85 |  26.38 |
  ---------+--------+--------+
         1 |    987 |    293 |   1280
           |  59.42 |  17.64 |  77.06
           |  77.11 |  22.89 |
           |  78.15 |  73.62 |
  ---------+--------+--------+
  Total        1263      398     1661
              76.04    23.96   100.00

         Frequency Missing = 29

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      3.5120    0.0609
 Likelihood Ratio Chi-Square    1      3.4351    0.0638
 Continuity Adj. Chi-Square     1      3.2605    0.0710
 Mantel-Haenszel Chi-Square     1      3.5099    0.0610
 Phi Coefficient                      -0.0460
 Contingency Coefficient               0.0459
 Cramer's V                           -0.0460


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)       276
           Left-sided Pr <= F          0.0365
           Right-sided Pr >= F         0.9730

           Table Probability (P)       0.0095
           Two-sided Pr <= P           0.0650

              Effective Sample Size = 1661

*/
