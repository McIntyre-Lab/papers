
/* Get counts and test for T1D enrichment */


libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname con '/home/jrbnewman/concannon/sas_data';



/* get counts for FDR 5% */

data eqtl_summary;
  set eqtl.eqtl_results_summary;
  keep snp_id feature_id feature_type flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05;
  if flag_cd4_fdr05=. then flag_cd4_fdr05=0;
    if flag_cd8_fdr05=. then flag_cd8_fdr05=0;
      if flag_cd19_fdr05=. then flag_cd19_fdr05=0;
  run;

proc sort data=eqtl_summary nodup;
    by snp_id feature_id;
    run;
        
proc freq data=eqtl_summary noprint;
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_all;
   run;
   
proc freq data=eqtl_summary noprint;
   where feature_type='exon';
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_exon;
   run;
   
proc freq data=eqtl_summary noprint;
  where feature_type='AS' or feature_type='IR';
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_splicing;
   run;
   
   proc freq data=eqtl_summary noprint;
   where feature_type='AS';
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_as;
   run;
   
   proc freq data=eqtl_summary noprint;
   where feature_type='IR';
   tables flag_cd4_fdr05*flag_cd8_fdr05*flag_cd19_fdr05 / out=eqtl_sig_fdr05_ir;
   run;
   
proc print data=eqtl_sig_fdr05_all;
run; quit;


/*
                          flag_
flag_cd4_    flag_cd8_    cd19_
  fdr05        fdr05      fdr05     COUNT    PERCENT

    0            0          0      435514    98.2512
    0            0          1        2063     0.4654
    0            1          0        1936     0.4368
    0            1          1         219     0.0494
    1            0          0        1911     0.4311
    1            0          1         326     0.0735
    1            1          0         333     0.0751
    1            1          1         964     0.2175



*/

proc print data=eqtl_sig_fdr05_exon;
run; quit;

/*
                          flag_
flag_cd4_    flag_cd8_    cd19_
  fdr05        fdr05      fdr05     COUNT    PERCENT

    0            0          0      189652    98.3448
    0            0          1         740     0.3837
    0            1          0         751     0.3894
    0            1          1         111     0.0576
    1            0          0         840     0.4356
    1            0          1         151     0.0783
    1            1          0         174     0.0902
    1            1          1         425     0.2204



*/

proc print data=eqtl_sig_fdr05_splicing;
run; quit;

/*
                            flag_
  flag_cd4_    flag_cd8_    cd19_
    fdr05        fdr05      fdr05     COUNT    PERCENT

      0            0          0      245862    98.1791
      0            0          1        1323     0.5283
      0            1          0        1185     0.4732
      0            1          1         108     0.0431
      1            0          0        1071     0.4277
      1            0          1         175     0.0699
      1            1          0         159     0.0635
      1            1          1         539     0.2152



*/

proc print data=eqtl_sig_fdr05_as;
run; quit;

/*

                          flag_
flag_cd4_    flag_cd8_    cd19_
  fdr05        fdr05      fdr05     COUNT    PERCENT

    0            0          0      202020    97.9976
    0            0          1        1183     0.5739
    0            1          0        1060     0.5142
    0            1          1         105     0.0509
    1            0          0         967     0.4691
    1            0          1         160     0.0776
    1            1          0         146     0.0708
    1            1          1         507     0.2459

*/

proc print data=eqtl_sig_fdr05_ir;
run; quit;

/*
                           flag_
 flag_cd4_    flag_cd8_    cd19_
   fdr05        fdr05      fdr05    COUNT    PERCENT

     0            0          0      43842    99.0243
     0            0          1        140     0.3162
     0            1          0        125     0.2823
     0            1          1          3     0.0068
     1            0          0        104     0.2349
     1            0          1         15     0.0339
     1            1          0         13     0.0294
     1            1          1         32     0.0723


*/

/* Enrichments for T1D genes */

data eqtl_summary_gene;
  set eqtl.eqtl_results_summary;
  if flag_cd4_fdr05=1 or flag_cd8_fdr05=1 or flag_cd19_fdr05=1 then flag_any_sig=1;
  else flag_any_sig=0;
    keep snp_id feature_id gene_id feature_type flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05 flag_any_sig;
  if flag_cd4_fdr05=. then flag_cd4_fdr05=0;
    if flag_cd8_fdr05=. then flag_cd8_fdr05=0;
      if flag_cd19_fdr05=. then flag_cd19_fdr05=0;
  run;

/* get gene info for fusions */

data fusion2gene;
   length fusion_id $2550.;
   format fusion_id $2475.;
   set fus.unique_info_fusions_si;
   rename fusion_id=feature_id;
   keep fusion_id gene_id;
run;

proc sort data=fusion2gene;
   by feature_id;
proc sort data=eqtl_summary_gene;
   by feature_id;
run;

data eqtl_summary_gene_fus eqtl_summary_gene_splice;
   set eqtl_summary_gene;
   if feature_type='exon' then output eqtl_summary_gene_fus;
   else output eqtl_summary_gene_splice;
run;



data eqtl_summary_gene_fus2;
   set eqtl_summary_gene_fus;
   drop gene_id;
run;


proc sort data=eqtl_summary_gene_fus2;
   by feature_id;
run;


data eqtl_summary_gene2 noeqtl oops;
   merge fusion2gene (in=in1) eqtl_summary_gene_fus2 (in=in2);
   by feature_id;
   if in1 and in2 then output eqtl_summary_gene2;
   else if in1 then output noeqtl;
   else output oops;
run;

data eqtl_summary_gene3;
   set eqtl_summary_gene2 eqtl_summary_gene_splice;
run;


data immunogenes;
   set con.immunogene_flags;
   keep gene_id flag_diabetes_gene;
   run;
   
proc sort data=immunogenes;
   by gene_id;
proc sort data=eqtl_summary_gene3 nodup;
   by gene_id;
   run;


data eqtl_summary_w_geneflags;
   merge eqtl_summary_gene3 (in=in1) immunogenes (in=in2);
   by gene_id;
   if in1 and in2 then output;
   else if in1 then do;
      flag_diabetes_gene=0; output; end;
   run;
   
/* Enrichments! */

* T1D genes in any sig events - all eqtl types;
proc freq data=eqtl_summary_w_geneflags;
   tables flag_any_sig*flag_diabetes_gene / chisq;
   run;


/*
Table of flag_any_sig by flag_diabetes_gene

    flag_any_sig     flag_diabetes_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 | 414241 |  21273 | 435514
             |  93.45 |   4.80 |  98.25
             |  95.12 |   4.88 |
             |  98.24 |  98.50 |
    ---------+--------+--------+
           1 |   7427 |    325 |   7752
             |   1.68 |   0.07 |   1.75
             |  95.81 |   4.19 |
             |   1.76 |   1.50 |
    ---------+--------+--------+
    Total      421668    21598   443266
                95.13     4.87   100.00

               The SAS System    07:45 Thursday
Statistics for Table of flag_any_sig by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      7.8712    0.0050
  Likelihood Ratio Chi-Square    1      8.2385    0.0041
  Continuity Adj. Chi-Square     1      7.7226    0.0055
  Mantel-Haenszel Chi-Square     1      7.8712    0.0050
  Phi Coefficient                      -0.0042
  Contingency Coefficient               0.0042
  Cramer's V                           -0.0042


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    414241
            Left-sided Pr <= F          0.0023
            Right-sided Pr >= F         0.9981

            Table Probability (P)       0.0004
            Two-sided Pr <= P           0.0048

                   Sample Size = 443266




*/

* T1D genes in any sig events - exon eqtl types;
proc freq data=eqtl_summary_w_geneflags;
   where feature_type='exon';
   tables flag_any_sig*flag_diabetes_gene / chisq;
   run;

/*
 
 Table of flag_any_sig by flag_diabetes_gene

     flag_any_sig     flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 | 181193 |   8459 | 189652
              |  93.96 |   4.39 |  98.34
              |  95.54 |   4.46 |
              |  98.34 |  98.47 |
     ---------+--------+--------+
            1 |   3061 |    131 |   3192
              |   1.59 |   0.07 |   1.66
              |  95.90 |   4.10 |
              |   1.66 |   1.53 |
     ---------+--------+--------+
     Total      184254     8590   192844
                 95.55     4.45   100.00

                The SAS System    07:45 Thursday,
 Statistics for Table of flag_any_sig by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      0.9362    0.3333
   Likelihood Ratio Chi-Square    1      0.9601    0.3272
   Continuity Adj. Chi-Square     1      0.8543    0.3553
   Mantel-Haenszel Chi-Square     1      0.9362    0.3333
   Phi Coefficient                      -0.0022
   Contingency Coefficient               0.0022
   Cramer's V                           -0.0022


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    181193
             Left-sided Pr <= F          0.1781
             Right-sided Pr >= F         0.8441

             Table Probability (P)       0.0222
             Two-sided Pr <= P           0.3634

                    Sample Size = 192844




*/

* T1D genes in any sig events - splicing eqtl types;
proc freq data=eqtl_summary_w_geneflags;
   where feature_type='AS' or feature_type='IR';
   tables flag_any_sig*flag_diabetes_gene / chisq;
   run;

/*
 Table of flag_any_sig by flag_diabetes_gene

     flag_any_sig     flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 | 233048 |  12814 | 245862
              |  93.06 |   5.12 |  98.18
              |  94.79 |   5.21 |
              |  98.16 |  98.51 |
     ---------+--------+--------+
            1 |   4366 |    194 |   4560
              |   1.74 |   0.08 |   1.82
              |  95.75 |   4.25 |
              |   1.84 |   1.49 |
     ---------+--------+--------+
     Total      237414    13008   250422
                 94.81     5.19   100.00

                The SAS System    07:45 Thursday, N


 Statistics for Table of flag_any_sig by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      8.3344    0.0039
   Likelihood Ratio Chi-Square    1      8.8483    0.0029
   Continuity Adj. Chi-Square     1      8.1411    0.0043
   Mantel-Haenszel Chi-Square     1      8.3343    0.0039
   Phi Coefficient                      -0.0058
   Contingency Coefficient               0.0058
   Cramer's V                           -0.0058


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    233048
             Left-sided Pr <= F          0.0017
             Right-sided Pr >= F         0.9987

             Table Probability (P)       0.0004
             Two-sided Pr <= P           0.0034

                    Sample Size = 250422


*/

* T1D genes in any sig events - as eqtl types;
proc freq data=eqtl_summary_w_geneflags;
   where feature_type='AS';
   tables flag_any_sig*flag_diabetes_gene / chisq;
   run;

/*
  Table of flag_any_sig by flag_diabetes_gene

      flag_any_sig     flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 | 192258 |   9762 | 202020
               |  93.26 |   4.74 |  98.00
               |  95.17 |   4.83 |
               |  97.99 |  98.19 |
      ---------+--------+--------+
             1 |   3948 |    180 |   4128
               |   1.92 |   0.09 |   2.00
               |  95.64 |   4.36 |
               |   2.01 |   1.81 |
      ---------+--------+--------+
      Total      196206     9942   206148
                  95.18     4.82   100.00

                 The SAS System    07:45 Thur

               The FREQ Procedure

Statistics for Table of flag_any_sig by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      1.9612    0.1614
  Likelihood Ratio Chi-Square    1      2.0224    0.1550
  Continuity Adj. Chi-Square     1      1.8597    0.1727
  Mantel-Haenszel Chi-Square     1      1.9612    0.1614
  Phi Coefficient                      -0.0031
  Contingency Coefficient               0.0031
  Cramer's V                           -0.0031


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    192258
            Left-sided Pr <= F          0.0848
            Right-sided Pr >= F         0.9264

            Table Probability (P)       0.0112
            Two-sided Pr <= P           0.1744

                   Sample Size = 206148




*/

* T1D genes in any sig events - ir eqtl types;
proc freq data=eqtl_summary_w_geneflags;
   where feature_type='IR';
   tables flag_any_sig*flag_diabetes_gene / chisq;
   run;

/*

 Table of flag_any_sig by flag_diabetes_gene

     flag_any_sig     flag_diabetes_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |  40790 |   3052 |  43842
              |  92.13 |   6.89 |  99.02
              |  93.04 |   6.96 |
              |  98.99 |  99.54 |
     ---------+--------+--------+
            1 |    418 |     14 |    432
              |   0.94 |   0.03 |   0.98
              |  96.76 |   3.24 |
              |   1.01 |   0.46 |
     ---------+--------+--------+
     Total       41208     3066    44274
                 93.07     6.93   100.00

                The SAS System    07:45 Thursday,

              The FREQ Procedure\

Statistics for Table of flag_any_sig by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      9.1876    0.0024
  Likelihood Ratio Chi-Square    1     11.2822    0.0008
  Continuity Adj. Chi-Square     1      8.6194    0.0033
  Mantel-Haenszel Chi-Square     1      9.1873    0.0024
  Phi Coefficient                      -0.0144
  Contingency Coefficient               0.0144
  Cramer's V                           -0.0144


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     40790
            Left-sided Pr <= F          0.0007
            Right-sided Pr >= F         0.9997

            Table Probability (P)       0.0004
            Two-sided Pr <= P           0.0015

                   Sample Size = 44274




*/

* T1D genes in any sig events - CD4 eqtls;
proc freq data=eqtl_summary_w_geneflags;
   tables flag_cd4_fdr05*flag_diabetes_gene / chisq;
   run;

/*

   Table of flag_cd4_fdr05 by flag_diabetes_gene

        flag_cd4_fdr05
                  flag_diabetes_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 | 418238 |  21494 | 439732
                 |  94.35 |   4.85 |  99.20
                 |  95.11 |   4.89 |
                 |  99.19 |  99.52 |
        ---------+--------+--------+
               1 |   3430 |    104 |   3534
                 |   0.77 |   0.02 |   0.80
                 |  97.06 |   2.94 |
                 |   0.81 |   0.48 |
        ---------+--------+--------+
        Total      421668    21598   443266
                    95.13     4.87   100.00

                   The SAS System    07:45 Thursday,

                 The FREQ Procedure
 Statistics for Table of flag_cd4_fdr05 by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     28.6177    <.0001
    Likelihood Ratio Chi-Square    1     33.1092    <.0001
    Continuity Adj. Chi-Square     1     28.1996    <.0001
    Mantel-Haenszel Chi-Square     1     28.6176    <.0001
    Phi Coefficient                      -0.0080
    Contingency Coefficient               0.0080
    Cramer's V                           -0.0080


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)    418238
              Left-sided Pr <= F          <.0001
              Right-sided Pr >= F         1.0000

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 443266


*/

* T1D genes in any sig events - CD8 eqtls;
proc freq data=eqtl_summary_w_geneflags;
   tables flag_cd8_fdr05*flag_diabetes_gene / chisq;
   run;

/*

 Table of flag_cd8_fdr05 by flag_diabetes_gene

      flag_cd8_fdr05
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 | 418372 |  21442 | 439814
               |  94.38 |   4.84 |  99.22
               |  95.12 |   4.88 |
               |  99.22 |  99.28 |
      ---------+--------+--------+
             1 |   3296 |    156 |   3452
               |   0.74 |   0.04 |   0.78
               |  95.48 |   4.52 |
               |   0.78 |   0.72 |
      ---------+--------+--------+
      Total      421668    21598   443266
                  95.13     4.87   100.00

                 The SAS System    07:45 Thursday,

 Statistics for Table of flag_cd8_fdr05 by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      0.9372    0.3330
    Likelihood Ratio Chi-Square    1      0.9593    0.3274
    Continuity Adj. Chi-Square     1      0.8619    0.3532
    Mantel-Haenszel Chi-Square     1      0.9372    0.3330
    Phi Coefficient                      -0.0015
    Contingency Coefficient               0.0015
    Cramer's V                           -0.0015


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)    418372
              Left-sided Pr <= F          0.1770
              Right-sided Pr >= F         0.8433

              Table Probability (P)       0.0203
              Two-sided Pr <= P           0.3612

                     Sample Size = 443266





*/


* T1D genes in any sig events - CD19 eqtls;
proc freq data=eqtl_summary_w_geneflags;
   tables flag_cd19_fdr05*flag_diabetes_gene / chisq;
   run;


/*
 Table of flag_cd19_fdr05 by flag_diabetes_gene

      flag_cd19_fdr05
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 | 418228 |  21466 | 439694
               |  94.35 |   4.84 |  99.19
               |  95.12 |   4.88 |
               |  99.18 |  99.39 |
      ---------+--------+--------+
             1 |   3440 |    132 |   3572
               |   0.78 |   0.03 |   0.81
               |  96.30 |   3.70 |
               |   0.82 |   0.61 |
      ---------+--------+--------+
      Total      421668    21598   443266
                  95.13     4.87   100.00

                 The SAS System    07:45 Thursday, No

               The FREQ Procedure


Statistics for Table of flag_cd19_fdr05 by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     10.7639    0.0010
    Likelihood Ratio Chi-Square    1     11.6954    0.0006
    Continuity Adj. Chi-Square     1     10.5094    0.0012
    Mantel-Haenszel Chi-Square     1     10.7638    0.0010
    Phi Coefficient                      -0.0049
    Contingency Coefficient               0.0049
    Cramer's V                           -0.0049


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)    418228
              Left-sided Pr <= F          0.0004
              Right-sided Pr >= F         0.9997

              Table Probability (P)       0.0001
              Two-sided Pr <= P           0.0008

                     Sample Size = 443266





*/
   

/* get counts for FDR 20% */

data eqtl_summary_fdr20;
  set eqtl.eqtl_results_summary;
  if flag_cd19_fdr20 lt 0.2 then flag_cd19_fdr20=1; else flag_cd19_fdr20=0;
    if flag_cd4_fdr20 lt 0.2 then flag_cd4_fdr20=1; else flag_cd4_fdr20=0;
      if flag_cd8_fdr20 lt 0.2 then flag_cd8_fdr20=1; else flag_cd8_fdr20=0;
  keep snp_id feature_id feature_type flag_cd4_fdr20 flag_cd8_fdr20 flag_cd19_fdr20;
  run;

proc sort data=eqtl_summaryi_fdr20 nodup;
    by snp_id feature_id;
    run;
    
proc freq data=eqtl_summary_fdr20 noprint;
   tables flag_cd4_fdr20*flag_cd8_fdr20*flag_cd19_fdr20 / out=eqtl_sig_fdr20_all;
   run;
   
proc freq data=eqtl_summary_fdr20 noprint;
   where feature_type='exon';
   tables  flag_cd4_fdr20*flag_cd8_fdr20*flag_cd19_fdr20 / out=eqtl_sig_fdr20_exon;
   run;
   
proc freq data=eqtl_summary_fdr20 noprint;
   where feature_type='AS' or feature_type='IR';
   tables  flag_cd4_fdr20*flag_cd8_fdr20*flag_cd19_fdr20 / out=eqtl_sig_fdr20_splicing;
   run;
   
   proc freq data=eqtl_summary_fdr20 noprint;
   where feature_type='AS';
   tables flag_cd4_fdr20*flag_cd8_fdr20*flag_cd19_fdr20 / out=eqtl_sig_fdr20_as;
   run;
   
   proc freq data=eqtl_summary_fdr20 noprint;
   where feature_type='IR';
   tables flag_cd4_fdr20*flag_cd8_fdr20*flag_cd19_fdr20 / out=eqtl_sig_fdr20_ir;
   run;
   
proc print data=eqtl_sig_fdr20_all;
run; quit;

/*
                                  flag_
        flag_cd4_    flag_cd8_    cd19_
 Obs      fdr20        fdr20      fdr20     COUNT     PERCENT

  1         1            1          1      1193553      100


*/

proc print data=eqtl_sig_fdr20_exon;
run; quit;


/*

                                  flag_
        flag_cd4_    flag_cd8_    cd19_
 Obs      fdr20        fdr20      fdr20     COUNT    PERCENT

  1         1            1          1      536036      100




*/
proc print data=eqtl_sig_fdr20_splicing;
run; quit;

/*
                                  flag_
        flag_cd4_    flag_cd8_    cd19_
 Obs      fdr20        fdr20      fdr20     COUNT    PERCENT

  1         1            1          1      657517      100



*/

proc print data=eqtl_sig_fdr20_as;
run; quit;

/*

                                 flag_
       flag_cd4_    flag_cd8_    cd19_
Obs      fdr20        fdr20      fdr20     COUNT    PERCENT

 1         1            1          1      553258      100


*/

proc print data=eqtl_sig_fdr20_ir;
run; quit;


/*
                                  flag_
        flag_cd4_    flag_cd8_    cd19_
 Obs      fdr20        fdr20      fdr20     COUNT    PERCENT

  1         1            1          1      104259      100

*/


/* Get counts on gene - FDR 5% */

data eqtl_summary_by_gene_fus eqtl_summary_by_gene_splice; 
  set eqtl.eqtl_results_summary;

  if flag_cd4_fdr05=. then flag_cd4_fdr05=0;
  if flag_cd8_fdr05=. then flag_cd8_fdr05=0;
  if flag_cd19_fdr05=. then flag_cd19_fdr05=0;
      
  if flag_cd19_fdr20 lt 0.2 then flag_cd19_fdr20=1; else flag_cd19_fdr20=0;
  if flag_cd4_fdr20 lt 0.2 then flag_cd4_fdr20=1; else flag_cd4_fdr20=0;
  if flag_cd8_fdr20 lt 0.2 then flag_cd8_fdr20=1; else flag_cd8_fdr20=0;
  
  keep gene_id feature_id snp_id feature_type flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05 flag_cd4_fdr20 flag_cd8_fdr20 flag_cd19_fdr20;

  if feature_type='exon' then output eqtl_summary_by_gene_fus;
  else output eqtl_summary_by_gene_splice;
    run;


data fusion2gene;
   set fus.unique_info_fusions_si;
   rename fusion_id=feature_id;
run;

data eqtl_summary_by_gene_fus2;
   set eqtl_summary_by_gene_fus;
   drop gene_id;
run;

proc sort data=fusion2gene;
   by feature_id;
proc sort data=eqtl_summary_by_gene_fus2;
   by feature_id;
run;


data eqtl_summary_by_gene_fus3;
   merge eqtl_summary_by_gene_fus2 (in=in1) fusion2gene (in=in2);
   by feature_id;
   if in1;
run;


data eqtl_summary_by_gene;
   set eqtl_summary_by_gene_fus3 eqtl_summary_by_gene_splice;
run;

proc sort data=eqtl_summary_by_gene out=eqtl_summary_by_gene2 nodup;
   by gene_id feature_id snp_id ;
   run;

proc means data=eqtl_summary_by_gene2 noprint;
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05 flag_cd4_fdr20 flag_cd8_fdr20 flag_cd19_fdr20;
   output out=eqtls_summed_by_gene_all sum=;
run;

proc means data=eqtl_summary_by_gene2 noprint;
   where feature_type='exon';
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05 flag_cd4_fdr20 flag_cd8_fdr20 flag_cd19_fdr20;
   output out=eqtls_summed_by_gene_exon sum=;
run;

proc means data=eqtl_summary_by_gene2 noprint;
   where feature_type='AS' or feature_type='IR';
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05 flag_cd4_fdr20 flag_cd8_fdr20 flag_cd19_fdr20;
   output out=eqtls_summed_by_gene_splicing sum=;
run;

proc means data=eqtl_summary_by_gene2 noprint;
   where feature_type='AS';
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05 flag_cd4_fdr20 flag_cd8_fdr20 flag_cd19_fdr20;
   output out=eqtls_summed_by_gene_as sum=;
run;

proc means data=eqtl_summary_by_gene2 noprint;
   where feature_type='IR';
   by gene_id;
   var flag_cd4_fdr05 flag_cd8_fdr05 flag_cd19_fdr05 flag_cd4_fdr20 flag_cd8_fdr20 flag_cd19_fdr20;
   output out=eqtls_summed_by_gene_ir sum=;
run;


data eqtls_summed_by_gene_all2;
   set eqtls_summed_by_gene_all;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr20 gt 0 then flag_cd4_fdr20=1;
   if flag_cd8_fdr20 gt 0 then flag_cd8_fdr20=1;
   if flag_cd19_fdr20 gt 0 then flag_cd19_fdr20=1;
      
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;

data eqtls_summed_by_gene_exon2;
   set eqtls_summed_by_gene_exon;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr20 gt 0 then flag_cd4_fdr20=1;
   if flag_cd8_fdr20 gt 0 then flag_cd8_fdr20=1;
   if flag_cd19_fdr20 gt 0 then flag_cd19_fdr20=1;
      
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;


data eqtls_summed_by_gene_splicing2;
   set eqtls_summed_by_gene_splicing;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr20 gt 0 then flag_cd4_fdr20=1;
   if flag_cd8_fdr20 gt 0 then flag_cd8_fdr20=1;
   if flag_cd19_fdr20 gt 0 then flag_cd19_fdr20=1;
      
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;

data eqtls_summed_by_gene_as2;
   set eqtls_summed_by_gene_as;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr20 gt 0 then flag_cd4_fdr20=1;
   if flag_cd8_fdr20 gt 0 then flag_cd8_fdr20=1;
   if flag_cd19_fdr20 gt 0 then flag_cd19_fdr20=1;
   
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;

data eqtls_summed_by_gene_ir2;
   set eqtls_summed_by_gene_ir;
   if flag_cd4_fdr05 gt 0 then flag_cd4_fdr05=1;
   if flag_cd8_fdr05 gt 0 then flag_cd8_fdr05=1;
   if flag_cd19_fdr05 gt 0 then flag_cd19_fdr05=1;
   if flag_cd4_fdr20 gt 0 then flag_cd4_fdr20=1;
   if flag_cd8_fdr20 gt 0 then flag_cd8_fdr20=1;
   if flag_cd19_fdr20 gt 0 then flag_cd19_fdr20=1;
   
   if flag_cd4_fdr05=1 and flag_cd8_fdr05=1 and flag_cd19_fdr05=1 then flag_any_sig=1;
   else flag_any_sig=0;
run;


/* Merge in T1D gene flags */

data immunoflags;
   set con.immunogene_flags;
   keep gene_id flag_diabetes_gene;
   run;

proc sort data=eqtls_summed_by_gene_all2;
   by gene_id;
proc sort data=eqtls_summed_by_gene_exon2;
   by gene_id;
proc sort data=eqtls_summed_by_gene_splicing2;
   by gene_id;
proc sort data=eqtls_summed_by_gene_as2;
   by gene_id;
proc sort data=eqtls_summed_by_gene_ir2;
   by gene_id;
proc sort data=immunoflags;
   by gene_id;
   run;



data eqtls_by_gene_all;
   merge eqtls_summed_by_gene_all2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;
   
data eqtls_by_gene_exon;
   merge eqtls_summed_by_gene_exon2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;
   
   
data eqtls_by_gene_splicing;
   merge eqtls_summed_by_gene_splicing2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;

data eqtls_by_gene_as;
   merge eqtls_summed_by_gene_as2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;

data eqtls_by_gene_ir;
   merge eqtls_summed_by_gene_ir2 (in=in1) immunoflags;
   by gene_id;
   if in1;
   run;
   
   
/* Get counts for each subset */

proc freq data=eqtls_by_gene_all noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_all;
   run;

proc freq data=eqtls_by_gene_exon noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_exon;
   run;
   
proc freq data=eqtls_by_gene_splicing noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_splicing;
   run;

proc freq data=eqtls_by_gene_as noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_as;
   run;

proc freq data=eqtls_by_gene_ir noprint;
   tables flag_cd19_fdr05*flag_cd4_fdr05*flag_cd8_fdr05 / out=eqtl_counts_by_gene_ir;
   run;



proc freq data=eqtls_by_gene_all noprint;
   tables flag_cd19_fdr20*flag_cd4_fdr20*flag_cd8_fdr20 / out=eqtl_counts_by_gene_all20;
   run;

proc freq data=eqtls_by_gene_exon noprint;
   tables flag_cd19_fdr20*flag_cd4_fdr20*flag_cd8_fdr20 / out=eqtl_counts_by_gene_exon20;
   run;
   
proc freq data=eqtls_by_gene_splicing noprint;
   tables flag_cd19_fdr20*flag_cd4_fdr20*flag_cd8_fdr20 / out=eqtl_counts_by_gene_splicing20;
   run;

proc freq data=eqtls_by_gene_as noprint;
   tables flag_cd19_fdr20*flag_cd4_fdr20*flag_cd8_fdr20 / out=eqtl_counts_by_gene_as20;
   run;

proc freq data=eqtls_by_gene_ir noprint;
   tables flag_cd19_fdr20*flag_cd4_fdr20*flag_cd8_fdr20 / out=eqtl_counts_by_gene_ir20;
   run;


proc print data=eqtl_counts_by_gene_all; run; quit;

/* 
        flag_
        cd19_    flag_cd4_    flag_cd8_
 Obs    fdr05      fdr05        fdr05      COUNT    PERCENT

  1       0          0            0         1161    65.6674
  2       0          0            1           89     5.0339
  3       0          1            0           91     5.1471
  4       0          1            1           56     3.1674
  5       1          0            0          101     5.7127
  6       1          0            1           46     2.6018
  7       1          1            0           52     2.9412
  8       1          1            1          172     9.7285




*/

proc print data=eqtl_counts_by_gene_exon; run; quit;


/* 
        flag_
        cd19_    flag_cd4_    flag_cd8_
 Obs    fdr05      fdr05        fdr05      COUNT    PERCENT

  1       0          0            0         1317    75.6462
  2       0          0            1           82     4.7099
  3       0          1            0           71     4.0781
  4       0          1            1           43     2.4698
  5       1          0            0           82     4.7099
  6       1          0            1           22     1.2636
  7       1          1            0           33     1.8955
  8       1          1            1           91     5.2269




*/

proc print data=eqtl_counts_by_gene_splicing; run; quit;
/* 
        flag_
        cd19_    flag_cd4_    flag_cd8_
 Obs    fdr05      fdr05        fdr05      COUNT    PERCENT

  1       0          0            0         585     57.8063
  2       0          0            1          63      6.2253
  3       0          1            0          59      5.8300
  4       0          1            1          39      3.8538
  5       1          0            0          73      7.2134
  6       1          0            1          39      3.8538
  7       1          1            0          44      4.3478
  8       1          1            1         110     10.8696




*/

proc print data=eqtl_counts_by_gene_as; run; quit;
/* 
        flag_
        cd19_    flag_cd4_    flag_cd8_
 Obs    fdr05      fdr05        fdr05      COUNT    PERCENT

  1       0          0            0         575     59.3395
  2       0          0            1          55      5.6760
  3       0          1            0          54      5.5728
  4       0          1            1          39      4.0248
  5       1          0            0          71      7.3271
  6       1          0            1          35      3.6120
  7       1          1            0          43      4.4376
  8       1          1            1          97     10.0103



*/

proc print data=eqtl_counts_by_gene_ir; run; quit;
/* 
        flag_
        cd19_    flag_cd4_    flag_cd8_
 Obs    fdr05      fdr05        fdr05      COUNT    PERCENT

  1       0          0            0         681     84.0741
  2       0          0            1          38      4.6914
  3       0          1            0          27      3.3333
  4       0          1            1           7      0.8642
  5       1          0            0          30      3.7037
  6       1          0            1           4      0.4938
  7       1          1            0           9      1.1111
  8       1          1            1          14      1.7284


*/


