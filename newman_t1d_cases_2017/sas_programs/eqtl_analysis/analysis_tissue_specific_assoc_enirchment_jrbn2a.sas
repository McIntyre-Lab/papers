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
able of flag_tissue_spec_assoc by flag_exon_eqtl

      flag_tissue_spec_assoc
                flag_exon_eqtl

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |    940 |    870 |   1810
               |  13.11 |  12.14 |  25.25
               |  51.93 |  48.07 |
               |  22.96 |  28.30 |
      ---------+--------+--------+
             1 |   3154 |   2204 |   5358
               |  44.00 |  30.75 |  74.75
               |  58.87 |  41.13 |
               |  77.04 |  71.70 |
      ---------+--------+--------+
      Total        4094     3074     7168
                  57.11    42.89   100.00

Statistics for Table of flag_tissue_spec_assoc by flag_exon_eqtl

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     26.5392    <.0001
     Likelihood Ratio Chi-Square    1     26.4025    <.0001
     Continuity Adj. Chi-Square     1     26.2570    <.0001
     Mantel-Haenszel Chi-Square     1     26.5355    <.0001
     Phi Coefficient                      -0.0608
     Contingency Coefficient               0.0607
     Cramer's V                           -0.0608

 Statistics for Table of flag_tissue_spec_assoc by flag_exon_eqtl

                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)       940
                Left-sided Pr <= F          <.0001
                Right-sided Pr >= F         1.0000

                Table Probability (P)       <.0001
                Two-sided Pr <= P           <.0001

                   Effective Sample Size = 7168
                    Frequency Missing = 349501

              WARNING: 98% of the data are missing.


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
              0 |    931 |    879 |   1810
                |  12.99 |  12.26 |  25.25
                |  51.44 |  48.56 |
                |  27.12 |  23.53 |
       ---------+--------+--------+
              1 |   2502 |   2856 |   5358
                |  34.91 |  39.84 |  74.75
                |  46.70 |  53.30 |
                |  72.88 |  76.47 |
       ---------+--------+--------+
       Total        3433     3735     7168
                   47.89    52.11   100.00

            Frequency Missing = 349501


 Statistics for Table of flag_tissue_spec_assoc by flag_junc_eqtl

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1     12.1803    0.0005
      Likelihood Ratio Chi-Square    1     12.1722    0.0005
      Continuity Adj. Chi-Square     1     11.9912    0.0005
      Mantel-Haenszel Chi-Square     1     12.1786    0.0005
      Phi Coefficient                       0.0412
      Contingency Coefficient               0.0412
      Cramer's V                            0.0412

 Statistics for Table of flag_tissue_spec_assoc by flag_junc_eqtl

                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)       931
                Left-sided Pr <= F          0.9998
                Right-sided Pr >= F         0.0003

                Table Probability (P)       <.0001
                Two-sided Pr <= P           0.0005

                   Effective Sample Size = 7168
                    Frequency Missing = 349501

              WARNING: 98% of the data are missing.


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
                0 |   1749 |     61 |   1810
                  |  24.40 |   0.85 |  25.25
                  |  96.63 |   3.37 |
                  |  25.69 |  16.99 |
         ---------+--------+--------+
                1 |   5060 |    298 |   5358
                  |  70.59 |   4.16 |  74.75
                  |  94.44 |   5.56 |
                  |  74.31 |  83.01 |
         ---------+--------+--------+
         Total        6809      359     7168
                     94.99     5.01   100.00

              Frequency Missing = 349501

                    The SAS System       10:41 Thursday,

 Statistics for Table of flag_tissue_spec_assoc by flag_ir_eqtl

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1     13.6593    0.0002
     Likelihood Ratio Chi-Square    1     14.8164    0.0001
     Continuity Adj. Chi-Square     1     13.2025    0.0003
     Mantel-Haenszel Chi-Square     1     13.6574    0.0002
     Phi Coefficient                       0.0437
     Contingency Coefficient               0.0436
     Cramer's V                            0.0437

 Statistics for Table of flag_tissue_spec_assoc by flag_ir_eqtl

                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)      1749
               Left-sided Pr <= F          1.0000
               Right-sided Pr >= F         <.0001

               Table Probability (P)       <.0001
               Two-sided Pr <= P           0.0001

                  Effective Sample Size = 7168
                   Frequency Missing = 349501

             WARNING: 98% of the data are missing.


*/
