/* Set libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';


/* Assoc with tissue-specific features */

data cd4_spec cd8_spec cd19_spec;
  set eqtl.eqtl_results_summary_table;
  if (CD4_FDR_P ne . and CD4_FDR_P lt 0.05) and CD8_FDR_P eq . and CD19_FDR_P eq  . then output cd4_spec;
  if CD4_FDR_P eq . and (CD8_FDR_P ne . and CD8_FDR_P lt 0.05) and CD19_FDR_P eq . then output cd8_spec;
  if CD4_FDR_P eq . and CD8_FDR_P eq . and (CD19_FDR_P ne . and CD19_FDR_P lt 0.05) then output cd19_spec;
run;


/* Enrichment for IR-specificity */

data tissue_specific_enrichment;
   set eqtl.eqtl_results_summary_table;
   if CD4_FDR_P eq . or CD8_FDR_P eq . or CD19_FDR_P eq . then delete; *drop features with tissue-specific detection;
   if feature_type='IR' then flag_ir_eqtl=1;
   else flag_ir_eqtl=0;

   if feature_type='exon' then flag_exon_eqtl=1;
   else flag_exon_eqtl=0;

   if feature_type='Junc' then flag_junc_eqtl=1;
   else flag_junc_eqtl=0;

   if CD4_FDR_P ge 0.05 and CD8_FDR_P ge 0.05 and CD19_FDR_P ge 0.05 then flag_tissue_spec_assoc=.;
   else if CD4_FDR_P < 0.05 and CD8_FDR_P ge 0.05 and CD19_FDR_P ge 0.05 then flag_tissue_spec_assoc=1;
   else if CD4_FDR_P ge 0.05 and CD8_FDR_P < 0.05 and CD19_FDR_P ge 0.05 then flag_tissue_spec_assoc=1;
   else if CD4_FDR_P ge 0.05 and CD8_FDR_P ge 0.05 and CD19_FDR_P < 0.05 then flag_tissue_spec_assoc=1;
   else flag_tissue_spec_assoc=0;
run;


/* Exon-eqtls vs other: tissue-specificity of associations of features detected in all tissues */

proc freq data=tissue_specific_enrichment;
   tables flag_tissue_spec_assoc*flag_exon_eqtl / chisq;
run;

/*
Table of flag_tissue_spec_assoc by flag_exon_eqtl

       flag_tissue_spec_assoc
                 flag_exon_eqtl

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |    937 |    844 |   1781
                |  13.21 |  11.90 |  25.11
                |  52.61 |  47.39 |
                |  22.93 |  28.06 |
       ---------+--------+--------+
              1 |   3149 |   2164 |   5313
                |  44.39 |  30.50 |  74.89
                |  59.27 |  40.73 |
                |  77.07 |  71.94 |
       ---------+--------+--------+
       Total        4086     3008     7094
                   57.60    42.40   100.00

Statistics for Table of flag_tissue_spec_assoc by flag_exon_eqtl

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     24.2166    <.0001
     Likelihood Ratio Chi-Square    1     24.0840    <.0001
     Continuity Adj. Chi-Square     1     23.9447    <.0001
     Mantel-Haenszel Chi-Square     1     24.2132    <.0001
     Phi Coefficient                      -0.0584
     Contingency Coefficient               0.0583
     Cramer's V                           -0.0584

         Fisher's Exact Test
  ----------------------------------
  Cell (1,1) Frequency (F)       937
  Left-sided Pr <= F          <.0001
  Right-sided Pr >= F         1.0000

  Table Probability (P)       <.0001
  Two-sided Pr <= P           <.0001

     Effective Sample Size = 7094
      Frequency Missing = 347150
*/

/* Junction-eqtls vs other: tissue-specificity of associations of features detected in all tissues */

proc freq data=tissue_specific_enrichment;
   tables flag_tissue_spec_assoc*flag_junc_eqtl / chisq;
run;

/*
Table of flag_tissue_spec_assoc by flag_junc_eqtl

       flag_tissue_spec_assoc
                 flag_junc_eqtl

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |    905 |    876 |   1781
                |  12.76 |  12.35 |  25.11
                |  50.81 |  49.19 |
                |  26.88 |  23.50 |
       ---------+--------+--------+
              1 |   2462 |   2851 |   5313
                |  34.71 |  40.19 |  74.89
                |  46.34 |  53.66 |
                |  73.12 |  76.50 |
       ---------+--------+--------+
       Total        3367     3727     7094
                   47.46    52.54   100.00

 Statistics for Table of flag_tissue_spec_assoc by flag_junc_eqtl

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1     10.7121    0.0011
      Likelihood Ratio Chi-Square    1     10.7022    0.0011
      Continuity Adj. Chi-Square     1     10.5334    0.0012
      Mantel-Haenszel Chi-Square     1     10.7106    0.0011
      Phi Coefficient                       0.0389
      Contingency Coefficient               0.0388
      Cramer's V                            0.0389


             Fisher's Exact Test
      ----------------------------------
      Cell (1,1) Frequency (F)       905
      Left-sided Pr <= F          0.9995
      Right-sided Pr >= F         0.0006

      Table Probability (P)       0.0001
      Two-sided Pr <= P           0.0011

         Effective Sample Size = 7094
*/

/* IR-eqtls vs other: tissue-specificity of associations of features detected in all tissues */

proc freq data=tissue_specific_enrichment;
   tables flag_tissue_spec_assoc*flag_ir_eqtl / chisq;
run;

/*
Table of flag_tissue_spec_assoc by flag_ir_eqtl

      flag_tissue_spec_assoc
                flag_ir_eqtl

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   1720 |     61 |   1781
               |  24.25 |   0.86 |  25.11
               |  96.57 |   3.43 |
               |  25.54 |  16.99 |
      ---------+--------+--------+
             1 |   5015 |    298 |   5313
               |  70.69 |   4.20 |  74.89
               |  94.39 |   5.61 |
               |  74.46 |  83.01 |
      ---------+--------+--------+
      Total        6735      359     7094
                  94.94     5.06   100.00

 Statistics for Table of flag_tissue_spec_assoc by flag_ir_eqtl

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     13.2405    0.0003
     Likelihood Ratio Chi-Square    1     14.3487    0.0002
     Continuity Adj. Chi-Square     1     12.7899    0.0003
     Mantel-Haenszel Chi-Square     1     13.2387    0.0003
     Phi Coefficient                       0.0432
     Contingency Coefficient               0.0432
     Cramer's V                            0.0432

Statistics for Table of flag_tissue_spec_assoc by flag_ir_eqtl

                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      1720
              Left-sided Pr <= F          0.9999
              Right-sided Pr >= F         0.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           0.0002

                 Effective Sample Size = 7094
*/
