/* Enrichments for splicing events */

/* We want:

1. Enrichment of T1D genes/autoimmune genes in detected splicing events
2. IR tissue specificity (All genes, T1D genes, autoimmune genes)
3. Unannotated junc tissue specificity (All genes, T1D genes, autoimmune genes)

*/

libname splicing '/mnt/data/splicing/';
libname con '/home/jrbnewman/concannon/sas_data';

data immunogenes;
   set con.immunogene_flags;
run;

proc sort data=immunogenes;
   by gene_id;
proc sort data=splicing.splicing_results_w_annot_fdr;
   by gene_id;
run;

data splicing_results_w_geneflags;
    merge splicing.splicing_results_w_annot_fdr (in=in1) immunogenes;
    by gene_id;
    if in1;
run;

data splicing_results_w_geneflags2;
    set splicing_results_w_geneflags;
    if flag_cd19_on=1 or flag_cd4_on=1 or flag_cd8_on=1 then flag_any_on=1;
    else flag_any_on=0;

    if flag_cd19_on=1 and flag_cd4_on=0 and flag_cd8_on=0 then flag_cd19_specific=1;
    else flag_cd19_specific=0;

    if flag_cd19_on=0 and flag_cd4_on=1 and flag_cd8_on=0 then flag_cd4_specific=1;
    else flag_cd4_specific=0;

    if flag_cd19_on=0 and flag_cd4_on=0 and flag_cd8_on=1 then flag_cd8_specific=1;
    else flag_cd8_specific=0;

    if flag_cd19_specific=1 or flag_cd4_specific=1 or flag_cd8_specific=1 then flag_cell_specific=1;
    else flag_cell_specific=0;

    if flag_immuno_gene=. then flag_immuno_gene=0;
    if flag_diabetes_gene=. then flag_diabetes_gene=0;
run;



/* Enrichments */

/* Detected events - Autoimmune gene enrichment */

proc freq data=splicing_results_w_geneflags2;
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
         0 |5574349 | 638159 |6212508
           |  86.76 |   9.93 |  96.70
           |  89.73 |  10.27 |
           |  96.74 |  96.28 |
  ---------+--------+--------+
         1 | 187567 |  24668 | 212235
           |   2.92 |   0.38 |   3.30
           |  88.38 |  11.62 |
           |   3.26 |   3.72 |
  ---------+--------+--------+
  Total     5761916   662827  6424743
              89.68    10.32   100.00

             The SAS System    07:45 Thursd
 Statistics for Table of flag_any_on by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1    404.7208    <.0001
  Likelihood Ratio Chi-Square    1    390.9953    <.0001
  Continuity Adj. Chi-Square     1    404.5748    <.0001
  Mantel-Haenszel Chi-Square     1    404.7207    <.0001
  Phi Coefficient                       0.0079
  Contingency Coefficient               0.0079
  Cramer's V                            0.0079


                   Fisher's Exact Test
           -----------------------------------
           Cell (1,1) Frequency (F)    5574349
           Left-sided Pr <= F           1.0000
           Right-sided Pr >= F          <.0001

           Table Probability (P)        <.0001
           Two-sided Pr <= P            <.0001

                  Sample Size = 6424743


*/

/* Detected events - T1D gene enrichment */

proc freq data=splicing_results_w_geneflags2;
     tables flag_any_on*flag_diabetes_gene  / chisq;
run;



/*
Table of flag_any_on by flag_diabetes_gene

   flag_any_on     flag_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |6204795 |   7713 |6212508
            |  96.58 |   0.12 |  96.70
            |  99.88 |   0.12 |
            |  96.70 |  94.33 |
   ---------+--------+--------+
          1 | 211771 |    464 | 212235
            |   3.30 |   0.01 |   3.30
            |  99.78 |   0.22 |
            |   3.30 |   5.67 |
   ---------+--------+--------+
   Total     6416566     8177  6424743
               99.87     0.13   100.00

              The SAS System    07:45 Thur
Statistics for Table of flag_any_on by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1    144.0976    <.0001
  Likelihood Ratio Chi-Square    1    119.2831    <.0001
  Continuity Adj. Chi-Square     1    143.3554    <.0001
  Mantel-Haenszel Chi-Square     1    144.0976    <.0001
  Phi Coefficient                       0.0047
  Contingency Coefficient               0.0047
  Cramer's V                            0.0047


                   Fisher's Exact Test
           -----------------------------------
           Cell (1,1) Frequency (F)    6204795
           Left-sided Pr <= F           1.0000
           Right-sided Pr >= F          <.0001

           Table Probability (P)        <.0001
           Two-sided Pr <= P            <.0001

                  Sample Size = 6424743

*/


/* IR event specificity - all genes */

proc freq data=splicing_results_w_geneflags2;
   where flag_any_on=1;
   tables flag_cell_specific*flag_intron_retention / chisq ;
run;


/*

 Table of flag_cell_specific by flag_intron_retention

         flag_cell_specific
                   flag_intron_retention

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 | 154948 |  36343 | 191291
                  |  73.01 |  17.12 |  90.13
                  |  81.00 |  19.00 |
                  |  91.59 |  84.41 |
         ---------+--------+--------+
                1 |  14232 |   6712 |  20944
                  |   6.71 |   3.16 |   9.87
                  |  67.95 |  32.05 |
                  |   8.41 |  15.59 |
         ---------+--------+--------+
         Total      169180    43055   212235
                     79.71    20.29   100.00

                    The SAS System    07:45 Thursday, N


 Statistics for Table of flag_cell_specific by flag_intron_retention

        Statistic                     DF       Value      Prob
        ------------------------------------------------------
        Chi-Square                     1   1987.5764    <.0001
        Likelihood Ratio Chi-Square    1   1793.9492    <.0001
        Continuity Adj. Chi-Square     1   1986.7696    <.0001
        Mantel-Haenszel Chi-Square     1   1987.5670    <.0001
        Phi Coefficient                       0.0968
        Contingency Coefficient               0.0963
        Cramer's V                            0.0968


                         Fisher's Exact Test
                  ----------------------------------
                  Cell (1,1) Frequency (F)    154948
                  Left-sided Pr <= F          1.0000
                  Right-sided Pr >= F         <.0001

                  Table Probability (P)       <.0001
                  Two-sided Pr <= P           <.0001

                         Sample Size = 212235



*/

/* Unannotated junction specificity - all genes */

proc freq data=splicing_results_w_geneflags2;
   where flag_any_on=1 and flag_intron_retention=0;
   tables flag_cell_specific*flag_junction_annotated / chisq ;
run;


/*

 Table of flag_cell_specific by flag_junction_annotated

          flag_cell_specific
                    flag_junction_annotated

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       0|       1|  Total
          ---------+--------+--------+
                 0 |  34101 | 120847 | 154948
                   |  20.16 |  71.43 |  91.59
                   |  22.01 |  77.99 |
                   |  88.02 |  92.65 |
          ---------+--------+--------+
                 1 |   4640 |   9592 |  14232
                   |   2.74 |   5.67 |   8.41
                   |  32.60 |  67.40 |
                   |  11.98 |   7.35 |
          ---------+--------+--------+
          Total       38741   130439   169180
                      22.90    77.10   100.00

                     The SAS System    07:45 Thursday, Novembe


 Statistics for Table of flag_cell_specific by flag_junction_annotated

         Statistic                     DF       Value      Prob
         ------------------------------------------------------
         Chi-Square                     1    828.6827    <.0001
         Likelihood Ratio Chi-Square    1    768.4947    <.0001
         Continuity Adj. Chi-Square     1    828.0827    <.0001
         Mantel-Haenszel Chi-Square     1    828.6778    <.0001
         Phi Coefficient                      -0.0700
         Contingency Coefficient               0.0698
         Cramer's V                           -0.0700


                          Fisher's Exact Test
                   ----------------------------------
                   Cell (1,1) Frequency (F)     34101
                   Left-sided Pr <= F          <.0001
                   Right-sided Pr >= F         1.0000

                   Table Probability (P)       <.0001
                   Two-sided Pr <= P           <.0001

                          Sample Size = 169180

*/




/* IR event specificity - autoimmune genes */

proc freq data=splicing_results_w_geneflags2;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_cell_specific*flag_intron_retention / chisq ;
run;


/*
 Table of flag_cell_specific by flag_intron_retention

         flag_cell_specific
                   flag_intron_retention

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |  18406 |   3871 |  22277
                  |  74.61 |  15.69 |  90.31
                  |  82.62 |  17.38 |
                  |  92.26 |  82.05 |
         ---------+--------+--------+
                1 |   1544 |    847 |   2391
                  |   6.26 |   3.43 |   9.69
                  |  64.58 |  35.42 |
                  |   7.74 |  17.95 |
         ---------+--------+--------+
         Total       19950     4718    24668
                     80.87    19.13   100.00

                    The SAS System    07:45 Thursday, No

                  The FREQ Procedure

Statistics for Table of flag_cell_specific by flag_intron_retention

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1    454.6948    <.0001
       Likelihood Ratio Chi-Square    1    394.2802    <.0001
       Continuity Adj. Chi-Square     1    453.5287    <.0001
       Mantel-Haenszel Chi-Square     1    454.6763    <.0001
       Phi Coefficient                       0.1358
       Contingency Coefficient               0.1345
       Cramer's V                            0.1358


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)     18406
                 Left-sided Pr <= F          1.0000
                 Right-sided Pr >= F         <.0001

                 Table Probability (P)       <.0001
                 Two-sided Pr <= P           <.0001

                        Sample Size = 24668






*/

/* Unannotated junction specificity - autoimmune genes */

proc freq data=splicing_results_w_geneflags2;
   where flag_any_on=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_cell_specific*flag_junction_annotated / chisq ;
run;


/*

Table of flag_cell_specific by flag_junction_annotated

         flag_cell_specific
                   flag_junction_annotated

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |   3937 |  14469 |  18406
                  |  19.73 |  72.53 |  92.26
                  |  21.39 |  78.61 |
                  |  88.37 |  93.38 |
         ---------+--------+--------+
                1 |    518 |   1026 |   1544
                  |   2.60 |   5.14 |   7.74
                  |  33.55 |  66.45 |
                  |  11.63 |   6.62 |
         ---------+--------+--------+
         Total        4455    15495    19950
                     22.33    77.67   100.00

                    The SAS System    07:45 Thursday, Novemb


Statistics for Table of flag_cell_specific by flag_junction_annotated

        Statistic                     DF       Value      Prob
        ------------------------------------------------------
        Chi-Square                     1    121.4336    <.0001
        Likelihood Ratio Chi-Square    1    111.0880    <.0001
        Continuity Adj. Chi-Square     1    120.7335    <.0001
        Mantel-Haenszel Chi-Square     1    121.4275    <.0001
        Phi Coefficient                      -0.0780
        Contingency Coefficient               0.0778
        Cramer's V                           -0.0780


                         Fisher's Exact Test
                  ----------------------------------
                  Cell (1,1) Frequency (F)      3937
                  Left-sided Pr <= F          <.0001
                  Right-sided Pr >= F         1.0000

                  Table Probability (P)       <.0001
                  Two-sided Pr <= P           <.0001

                         Sample Size = 19950


*/



/* IR event specificity - T1D genes */

proc freq data=splicing_results_w_geneflags2;
   where flag_any_on=1 and flag_diabetes_gene=1;
   tables flag_cell_specific*flag_intron_retention / chisq ;
run;


/*
 Table of flag_cell_specific by flag_intron_retention

         flag_cell_specific
                   flag_intron_retention

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |    291 |     79 |    370
                  |  62.72 |  17.03 |  79.74
                  |  78.65 |  21.35 |
                  |  82.67 |  70.54 |
         ---------+--------+--------+
                1 |     61 |     33 |     94
                  |  13.15 |   7.11 |  20.26
                  |  64.89 |  35.11 |
                  |  17.33 |  29.46 |
         ---------+--------+--------+
         Total         352      112      464
                     75.86    24.14   100.00

                    The SAS System    07:45 Thursday, Nov


Statistics for Table of flag_cell_specific by flag_intron_retention

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      7.7448    0.0054
       Likelihood Ratio Chi-Square    1      7.2841    0.0070
       Continuity Adj. Chi-Square     1      7.0118    0.0081
       Mantel-Haenszel Chi-Square     1      7.7281    0.0054
       Phi Coefficient                       0.1292
       Contingency Coefficient               0.1281
       Cramer's V                            0.1292


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)       291
                 Left-sided Pr <= F          0.9977
                 Right-sided Pr >= F         0.0049

                 Table Probability (P)       0.0026
                 Two-sided Pr <= P           0.0069

                         Sample Size = 464


*/

/* Unannotated junction specificity - T1D genes */

proc freq data=splicing_results_w_geneflags2;
   where flag_any_on=1 and flag_intron_retention=0 and flag_diabetes_gene=1;
   tables flag_cell_specific*flag_junction_annotated / chisq ;
run;


/*
 Table of flag_cell_specific by flag_junction_annotated

          flag_cell_specific
                    flag_junction_annotated

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       0|       1|  Total
          ---------+--------+--------+
                 0 |     64 |    227 |    291
                   |  18.18 |  64.49 |  82.67
                   |  21.99 |  78.01 |
                   |  74.42 |  85.34 |
          ---------+--------+--------+
                 1 |     22 |     39 |     61
                   |   6.25 |  11.08 |  17.33
                   |  36.07 |  63.93 |
                   |  25.58 |  14.66 |
          ---------+--------+--------+
          Total          86      266      352
                      24.43    75.57   100.00

                     The SAS System    07:45 Thursday, Novem

                   The FREQ Procedure


     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1      5.4091    0.0200
     Likelihood Ratio Chi-Square    1      5.0558    0.0245
     Continuity Adj. Chi-Square     1      4.6737    0.0306
     Mantel-Haenszel Chi-Square     1      5.3937    0.0202
     Phi Coefficient                      -0.1240
     Contingency Coefficient               0.1230
     Cramer's V                           -0.1240


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)        64
               Left-sided Pr <= F          0.0175
               Right-sided Pr >= F         0.9922

               Table Probability (P)       0.0096
               Two-sided Pr <= P           0.0320

                       Sample Size = 352



*/

