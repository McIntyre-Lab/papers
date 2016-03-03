
/* Get counts and test for T1D enrichment */


libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';



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

     0            0          0      440102    98.2430
     0            0          1        2080     0.4643
     0            1          0        1977     0.4413
     0            1          1         219     0.0489
     1            0          0        1929     0.4306
     1            0          1         337     0.0752
     1            1          0         359     0.0801
     1            1          1         970     0.2165

*/

proc print data=eqtl_sig_fdr05_exon;
run; quit;

/*
                          flag_
flag_cd4_    flag_cd8_    cd19_
  fdr05        fdr05      fdr05     COUNT    PERCENT

    0            0          0      194055    98.3249
    0            0          1         756     0.3831
    0            1          0         788     0.3993
    0            1          1         110     0.0557
    1            0          0         861     0.4363
    1            0          1         160     0.0811
    1            1          0         200     0.1013
    1            1          1         431     0.2184



*/

proc print data=eqtl_sig_fdr05_as;
run; quit;

/*
                          flag_
flag_cd4_    flag_cd8_    cd19_
  fdr05        fdr05      fdr05     COUNT    PERCENT

    0            0          0      202108    97.9941
    0            0          1        1184     0.5741
    0            1          0        1064     0.5159
    0            1          1         106     0.0514
    1            0          0         968     0.4693
    1            0          1         162     0.0785
    1            1          0         146     0.0708
    1            1          1         507     0.2458

*/

proc print data=eqtl_sig_fdr05_ir;
run; quit;

/*
                             flag_
 flag_cd4_    flag_cd8_    cd19_
   fdr05        fdr05      fdr05    COUNT    PERCENT

     0            0          0      43939    99.0353
     0            0          1        140     0.3155
     0            1          0        125     0.2817
     0            1          1          3     0.0068
     1            0          0        100     0.2254
     1            0          1         15     0.0338
     1            1          0         13     0.0293
     1            1          1         32     0.0721

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

* T1D genes in any sig events - alSVl eqtl types;
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
            0 | 418207 |  21895 | 440102
              |  93.36 |   4.89 |  98.24
              |  95.03 |   4.97 |
              |  98.23 |  98.48 |
     ---------+--------+--------+
            1 |   7534 |    337 |   7871
              |   1.68 |   0.08 |   1.76
              |  95.72 |   4.28 |
              |   1.77 |   1.52 |
     ---------+--------+--------+
     Total      425741    22232   447973
                 95.04     4.96   100.00

 Statistics for Table of flag_any_sig by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      7.8838    0.0050
   Likelihood Ratio Chi-Square    1      8.2448    0.0041
   Continuity Adj. Chi-Square     1      7.7374    0.0054
   Mantel-Haenszel Chi-Square     1      7.8837    0.0050
   Phi Coefficient                      -0.0042
   Contingency Coefficient               0.0042
   Cramer's V                           -0.0042


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    418207
             Left-sided Pr <= F          0.0023
             Right-sided Pr >= F         0.9981

             Table Probability (P)       0.0004
             Two-sided Pr <= P           0.0047

                    Sample Size = 447973


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
            0 | 184978 |   9077 | 194055
              |  93.73 |   4.60 |  98.32
              |  95.32 |   4.68 |
              |  98.32 |  98.46 |
     ---------+--------+--------+
            1 |   3164 |    142 |   3306
              |   1.60 |   0.07 |   1.68
              |  95.70 |   4.30 |
              |   1.68 |   1.54 |
     ---------+--------+--------+
     Total      188142     9219   197361
                 95.33     4.67   100.00

                The SAS System       10:41 Thursday, Fe


Statistics for Table of flag_any_sig by flag_diabetes_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      1.0670    0.3016
  Likelihood Ratio Chi-Square    1      1.0949    0.2954
  Continuity Adj. Chi-Square     1      0.9829    0.3215
  Mantel-Haenszel Chi-Square     1      1.0670    0.3016
  Phi Coefficient                      -0.0023
  Contingency Coefficient               0.0023
  Cramer's V                           -0.0023


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    184978
            Left-sided Pr <= F          0.1607
            Right-sided Pr >= F         0.8592

            Table Probability (P)       0.0199
            Two-sided Pr <= P           0.3185

                   Sample Size = 197361

*/


* T1D genes in any sig events - as eqtl types;
proc freq data=eqtl_summary_w_geneflags;
   where feature_type='AS';
   tables flag_any_sig*flag_diabetes_gene / chisq;
   run;

/*
                The FREQ Procedure

  Table of flag_any_sig by flag_diabetes_gene

      flag_any_sig     flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 | 192347 |   9761 | 202108
               |  93.26 |   4.73 |  97.99
               |  95.17 |   4.83 |
               |  97.98 |  98.18 |
      ---------+--------+--------+
             1 |   3956 |    181 |   4137
               |   1.92 |   0.09 |   2.01
               |  95.62 |   4.38 |
               |   2.02 |   1.82 |
      ---------+--------+--------+
      Total      196303     9942   206245
                  95.18     4.82   100.00


 Statistics for Table of flag_any_sig by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      1.8248    0.1767
   Likelihood Ratio Chi-Square    1      1.8796    0.1704
   Continuity Adj. Chi-Square     1      1.7271    0.1888
   Mantel-Haenszel Chi-Square     1      1.8248    0.1767
   Phi Coefficient                      -0.0030
   Contingency Coefficient               0.0030
   Cramer's V                           -0.0030


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    192347
             Left-sided Pr <= F          0.0930
             Right-sided Pr >= F         0.9189

             Table Probability (P)       0.0120
             Two-sided Pr <= P           0.1868

                    Sample Size = 206245


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
           0 |  40882 |   3057 |  43939
             |  92.15 |   6.89 |  99.04
             |  93.04 |   6.96 |
             |  99.00 |  99.54 |
    ---------+--------+--------+
           1 |    414 |     14 |    428
             |   0.93 |   0.03 |   0.96
             |  96.73 |   3.27 |
             |   1.00 |   0.46 |
    ---------+--------+--------+
    Total       41296     3071    44367
                93.08     6.92   100.00

               The SAS System       10:41 Thurs

 Statistics for Table of flag_any_sig by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      8.9404    0.0028
   Likelihood Ratio Chi-Square    1     10.9538    0.0009
   Continuity Adj. Chi-Square     1      8.3774    0.0038
   Mantel-Haenszel Chi-Square     1      8.9402    0.0028
   Phi Coefficient                      -0.0142
   Contingency Coefficient               0.0142
   Cramer's V                           -0.0142


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     40882
             Left-sided Pr <= F          0.0008
             Right-sided Pr >= F         0.9997

             Table Probability (P)       0.0005
             Two-sided Pr <= P           0.0015

                    Sample Size = 44367



*/

* T1D genes in any sig events - CD4 eqtls;
proc freq data=eqtl_summary_w_geneflags;
   tables flag_cd4_fdr05*flag_diabetes_gene / chisq;
   run;

/*

  flag_cd4_fdr05
            flag_diabetes_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 422260 |  22118 | 444378
           |  94.26 |   4.94 |  99.20
           |  95.02 |   4.98 |
           |  99.18 |  99.49 |
  ---------+--------+--------+
         1 |   3481 |    114 |   3595
           |   0.78 |   0.03 |   0.80
           |  96.83 |   3.17 |
           |   0.82 |   0.51 |
  ---------+--------+--------+
  Total      425741    22232   447973
              95.04     4.96   100.00

  Statistics for Table of flag_cd4_fdr05 by flag_diabetes_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     24.6673    <.0001
     Likelihood Ratio Chi-Square    1     28.1085    <.0001
     Continuity Adj. Chi-Square     1     24.2859    <.0001
     Mantel-Haenszel Chi-Square     1     24.6673    <.0001
     Phi Coefficient                      -0.0074
     Contingency Coefficient               0.0074
     Cramer's V                           -0.0074


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)    422260
               Left-sided Pr <= F          <.0001
               Right-sided Pr >= F         1.0000

               Table Probability (P)       <.0001
               Two-sided Pr <= P           <.0001

                      Sample Size = 447973



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
             0 | 422376 |  22072 | 444448
               |  94.29 |   4.93 |  99.21
               |  95.03 |   4.97 |
               |  99.21 |  99.28 |
      ---------+--------+--------+
             1 |   3365 |    160 |   3525
               |   0.75 |   0.04 |   0.79
               |  95.46 |   4.54 |
               |   0.79 |   0.72 |
      ---------+--------+--------+
      Total      425741    22232   447973
                  95.04     4.96   100.00

                 The SAS System       10:41 Thu

 Statistics for Table of flag_cd8_fdr05 by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      1.3529    0.2448
    Likelihood Ratio Chi-Square    1      1.3908    0.2383
    Continuity Adj. Chi-Square     1      1.2639    0.2609
    Mantel-Haenszel Chi-Square     1      1.3529    0.2448
    Phi Coefficient                      -0.0017
    Contingency Coefficient               0.0017
    Cramer's V                           -0.0017


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)    422376
              Left-sided Pr <= F          0.1297
              Right-sided Pr >= F         0.8864

              Table Probability (P)       0.0162
              Two-sided Pr <= P           0.2588

                     Sample Size = 447973

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
             0 | 422268 |  22099 | 444367
               |  94.26 |   4.93 |  99.20
               |  95.03 |   4.97 |
               |  99.18 |  99.40 |
      ---------+--------+--------+
             1 |   3473 |    133 |   3606
               |   0.78 |   0.03 |   0.80
               |  96.31 |   3.69 |
               |   0.82 |   0.60 |
      ---------+--------+--------+
      Total      425741    22232   447973
                  95.04     4.96   100.00

Statistics for Table of flag_cd19_fdr05 by flag_diabetes_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     12.5198    0.0004
    Likelihood Ratio Chi-Square    1     13.6812    0.0002
    Continuity Adj. Chi-Square     1     12.2488    0.0005
    Mantel-Haenszel Chi-Square     1     12.5197    0.0004
    Phi Coefficient                      -0.0053
    Contingency Coefficient               0.0053
    Cramer's V                           -0.0053


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)    422268
              Left-sided Pr <= F          0.0001
              Right-sided Pr >= F         0.9999

              Table Probability (P)       <.0001
              Two-sided Pr <= P           0.0003

                     Sample Size = 447973

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
  Obs      fdr20        fdr20      fdr20     COUNT    PERCENT

   1         1            1          1      447973      100


*/

proc print data=eqtl_sig_fdr20_exon;
run; quit;


/*
                                 flag_
       flag_cd4_    flag_cd8_    cd19_
Obs      fdr20        fdr20      fdr20     COUNT    PERCENT

 1         1            1          1      197361      100



*/

proc print data=eqtl_sig_fdr20_as;
run; quit;

/*
                                  flag_
        flag_cd4_    flag_cd8_    cd19_
 Obs      fdr20        fdr20      fdr20     COUNT    PERCENT

  1         1            1          1      206245      100


*/

proc print data=eqtl_sig_fdr20_ir;
run; quit;


/*
                                    flag_
          flag_cd4_    flag_cd8_    cd19_
   Obs      fdr20        fdr20      fdr20    COUNT    PERCENT

    1         1            1          1      44367      100
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
  fdr05      fdr05        fdr05      COUNT    PERCENT

    0          0            0         1279    66.7537
    0          0            1           94     4.9061
    0          1            0           95     4.9582
    0          1            1           59     3.0793
    1          0            0          111     5.7933
    1          0            1           48     2.5052
    1          1            0           54     2.8184
    1          1            1          176     9.1858

*/

proc print data=eqtl_counts_by_gene_exon; run; quit;


/* 
 flag_
 cd19_    flag_cd4_    flag_cd8_
 fdr05      fdr05        fdr05      COUNT    PERCENT

   0          0            0         1440    76.0296
   0          0            1           87     4.5935
   0          1            0           75     3.9599
   0          1            1           46     2.4287
   1          0            0           92     4.8574
   1          0            1           23     1.2144
   1          1            0           36     1.9007
   1          1            1           95     5.0158

*/


proc print data=eqtl_counts_by_gene_as; run; quit;
/* 
 flag_
cd19_    flag_cd4_    flag_cd8_
fdr05      fdr05        fdr05      COUNT    PERCENT

  0          0            0         578     59.4650
  0          0            1          55      5.6584
  0          1            0          54      5.5556
  0          1            1          38      3.9095
  1          0            0          70      7.2016
  1          0            1          36      3.7037
  1          1            0          43      4.4239
  1          1            1          98     10.0823


*/

proc print data=eqtl_counts_by_gene_ir; run; quit;
/* 
flag_
cd19_    flag_cd4_    flag_cd8_
fdr05      fdr05        fdr05      COUNT    PERCENT

  0          0            0         683     84.1133
  0          0            1          38      4.6798
  0          1            0          27      3.3251
  0          1            1           7      0.8621
  1          0            0          30      3.6946
  1          0            1           4      0.4926
  1          1            0           9      1.1084
  1          1            1          14      1.7241

*/


