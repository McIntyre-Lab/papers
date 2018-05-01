libname event "!MCLAB/event_analysis/sas_data";
libname eventloc "/mnt/store/event_sandbox/sas_data";
libname hg19 "!PATCON/useful_human_data/aceview_hg19/fusions/sas_data";
libname con "!PATCON/sas_data";

/* Testing if the MEEs in the case-only data are more frequently
   unique, common, or constitutive exons. I am going to use the "mean" MEE to test this first */

data mee;
  set eventloc.t1d_exons_mee_rank_by_cell_mean;
  if flag_mee_cd4=1 and flag_mee_cd8=1 then flag_mee_tcells=1; else flag_mee_tcells=0;
  if flag_mee_cd4=1 and flag_mee_cd8=1 and flag_mee_cd19=1 then flag_mee_all=1; else flag_mee_all=0;
  keep fusion_id flag_mee_cd19 flag_mee_cd4 flag_mee_cd8 flag_mee_tcells flag_mee_all;
run;

data ucx_flags;
  set hg19.hg19_si_fusions_unique_flagged;
  if flag_common=1 or flag_constitutive=1 then flag_common_constit=1; else flag_common_constit=0;
  keep fusion_id flag_multigene flag_alternative flag_common flag_constitutive flag_common_constit;
run;

proc sort data=mee nodup;
  by fusion_id;
proc sort data=ucx_flags nodup;
   by fusion_id;
run;

data mme_ucx_flags;
  merge mee (in=in1) ucx_flags (in=in2);
  by fusion_id;
  if in1 and in2;
run;

proc freq data=mme_ucx_flags;
   tables flag_mee_cd4*flag_mee_cd8*flag_mee_cd19*flag_alternative*flag_common*flag_Constitutive / out=fus_count;
proc print data=fus_count;
run;

/*
  flag_      flag_       flag_        flag_        flag_        flag_
 mee_cd4    mee_cd8    mee_cd19    Alternative    Common    Constitutive    COUNT

    0          0           0            0            0            1           469
    0          0           0            0            1            0           710
    0          0           0            1            0            0         38835
    0          0           1            0            0            1            11
    0          0           1            0            1            0            13
    0          0           1            1            0            0           167
    0          1           0            0            0            1             7
    0          1           0            0            1            0             3
    0          1           0            1            0            0            67
    0          1           1            0            0            1             8
    0          1           1            0            1            0             5
    0          1           1            1            0            0            47
    1          0           0            0            0            1             9
    1          0           0            0            1            0             6
    1          0           0            1            0            0            54
    1          0           1            0            0            1             7
    1          0           1            0            1            0             2
    1          0           1            1            0            0            59
    1          1           0            0            0            1            10
    1          1           0            0            1            0             9
    1          1           0            1            0            0           163
    1          1           1            0            0            1          2307
    1          1           1            0            1            0           510
    1          1           1            1            0            0          2674
*/

proc freq data=mme_ucx_flags;
   tables flag_mee_cd4*flag_common_constit / chisq;
run;

/*
  flag_mee_cd4
            flag_common_constit

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  39116 |   1226 |  40342
           |  84.75 |   2.66 |  87.41
           |  96.96 |   3.04 |
           |  92.99 |  30.00 |
  ---------+--------+--------+
         1 |   2950 |   2860 |   5810
           |   6.39 |   6.20 |  12.59
           |  50.77 |  49.23 |
           |   7.01 |  70.00 |
  ---------+--------+--------+
  Total       42066     4086    46152
              91.15     8.85   100.00


 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1  13425.2938    <.0001
 Likelihood Ratio Chi-Square    1   8577.3292    <.0001
 Continuity Adj. Chi-Square     1  13419.5709    <.0001
 Mantel-Haenszel Chi-Square     1  13425.0029    <.0001
 Phi Coefficient                       0.5393
 Contingency Coefficient               0.4747
 Cramer's V                            0.5393


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     39116
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                  Sample Size = 46152

*/

proc freq data=mme_ucx_flags;
   tables flag_mee_cd8*flag_common_constit / chisq;
run;

/*
  flag_mee_cd8
            flag_common_constit

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  39115 |   1227 |  40342
           |  84.75 |   2.66 |  87.41
           |  96.96 |   3.04 |
           |  92.98 |  30.03 |
  ---------+--------+--------+
         1 |   2951 |   2859 |   5810
           |   6.39 |   6.19 |  12.59
           |  50.79 |  49.21 |
           |   7.02 |  69.97 |
  ---------+--------+--------+
  Total       42066     4086    46152
              91.15     8.85   100.00

 Statistics for Table of flag_mee_cd8 by flag_common_constit

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1  13413.8491    <.0001
    Likelihood Ratio Chi-Square    1   8570.4672    <.0001
    Continuity Adj. Chi-Square     1  13408.1286    <.0001
    Mantel-Haenszel Chi-Square     1  13413.5585    <.0001
    Phi Coefficient                       0.5391
    Contingency Coefficient               0.4745
    Cramer's V                            0.5391


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     39115
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 46152
*/
proc freq data=mme_ucx_flags;
   tables flag_mee_cd19*flag_common_constit / chisq;
run;

/*
    flag_mee_cd19
              flag_common_constit

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  39119 |   1223 |  40342
             |  84.76 |   2.65 |  87.41
             |  96.97 |   3.03 |
             |  92.99 |  29.93 |
    ---------+--------+--------+
           1 |   2947 |   2863 |   5810
             |   6.39 |   6.20 |  12.59
             |  50.72 |  49.28 |
             |   7.01 |  70.07 |
    ---------+--------+--------+
    Total       42066     4086    46152
                91.15     8.85   100.00


  Statistics for Table of flag_mee_cd19 by flag_common_constit

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1  13459.6571    <.0001
     Likelihood Ratio Chi-Square    1   8597.9337    <.0001
     Continuity Adj. Chi-Square     1  13453.9269    <.0001
     Mantel-Haenszel Chi-Square     1  13459.3655    <.0001
     Phi Coefficient                       0.5400
     Contingency Coefficient               0.4752
     Cramer's V                            0.5400


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)     39119
               Left-sided Pr <= F          1.0000
               Right-sided Pr >= F         <.0001

               Table Probability (P)       <.0001
               Two-sided Pr <= P           <.0001

                      Sample Size = 46152

*/
proc freq data=mme_ucx_flags;
   tables flag_mee_tcells*flag_common_constit / chisq;
run;

/*
  flag_mee_tcells
            flag_common_constit

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  39229 |   1250 |  40479
           |  85.00 |   2.71 |  87.71
           |  96.91 |   3.09 |
           |  93.26 |  30.59 |
  ---------+--------+--------+
         1 |   2837 |   2836 |   5673
           |   6.15 |   6.14 |  12.29
           |  50.01 |  49.99 |
           |   6.74 |  69.41 |
  ---------+--------+--------+
  Total       42066     4086    46152
              91.15     8.85   100.00

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1  13564.6256    <.0001
   Likelihood Ratio Chi-Square    1   8591.5165    <.0001
   Continuity Adj. Chi-Square     1  13558.8138    <.0001
   Mantel-Haenszel Chi-Square     1  13564.3317    <.0001
   Phi Coefficient                       0.5421
   Contingency Coefficient               0.4766
   Cramer's V                            0.5421


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     39229
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 46152


*/
proc freq data=mme_ucx_flags;
   tables flag_mee_all*flag_common_constit / chisq;
run;

/*
  flag_mee_all
            flag_common_constit

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  39392 |   1269 |  40661
           |  85.35 |   2.75 |  88.10
           |  96.88 |   3.12 |
           |  93.64 |  31.06 |
  ---------+--------+--------+
         1 |   2674 |   2817 |   5491
           |   5.79 |   6.10 |  11.90
           |  48.70 |  51.30 |
           |   6.36 |  68.94 |
  ---------+--------+--------+
  Total       42066     4086    46152
              91.15     8.85   100.00

  Statistics for Table of flag_mee_all by flag_common_constit

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1  13917.0037    <.0001
     Likelihood Ratio Chi-Square    1   8705.3170    <.0001
     Continuity Adj. Chi-Square     1  13911.0336    <.0001
     Mantel-Haenszel Chi-Square     1  13916.7022    <.0001
     Phi Coefficient                       0.5491
     Contingency Coefficient               0.4813
     Cramer's V                            0.5491


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)     39392
               Left-sided Pr <= F          1.0000
               Right-sided Pr >= F         <.0001

               Table Probability (P)       <.0001
               Two-sided Pr <= P           <.0001

                      Sample Size = 46152
*/
proc freq data=mme_ucx_flags;
   tables flag_mee_cd4*flag_Constitutive / chisq;
run;

/*
   flag_mee_cd4     flag_Constitutive

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  39847 |    495 |  40342
            |  86.34 |   1.07 |  87.41
            |  98.77 |   1.23 |
            |  91.97 |  17.50 |
   ---------+--------+--------+
          1 |   3477 |   2333 |   5810
            |   7.53 |   5.06 |  12.59
            |  59.85 |  40.15 |
            |   8.03 |  82.50 |
   ---------+--------+--------+
   Total       43324     2828    46152
               93.87     6.13   100.00

 Statistics for Table of flag_mee_cd4 by flag_Constitutive

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1  13379.4420    <.0001
   Likelihood Ratio Chi-Square    1   8104.5859    <.0001
   Continuity Adj. Chi-Square     1  13372.6752    <.0001
   Mantel-Haenszel Chi-Square     1  13379.1521    <.0001
   Phi Coefficient                       0.5384
   Contingency Coefficient               0.4741
   Cramer's V                            0.5384


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     39847
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 46152



*/
proc freq data=mme_ucx_flags;
   tables flag_mee_cd8*flag_Constitutive / chisq;
run;

/*
  Table of flag_mee_cd8 by flag_Constitutive

     flag_mee_cd8     flag_Constitutive

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |  39846 |    496 |  40342
              |  86.34 |   1.07 |  87.41
              |  98.77 |   1.23 |
              |  91.97 |  17.54 |
     ---------+--------+--------+
            1 |   3478 |   2332 |   5810
              |   7.54 |   5.05 |  12.59
              |  59.86 |  40.14 |
              |   8.03 |  82.46 |
     ---------+--------+--------+
     Total       43324     2828    46152
                 93.87     6.13   100.00


  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1  13365.9102    <.0001
  Likelihood Ratio Chi-Square    1   8096.6102    <.0001
  Continuity Adj. Chi-Square     1  13359.1469    <.0001
  Mantel-Haenszel Chi-Square     1  13365.6206    <.0001
  Phi Coefficient                       0.5382
  Contingency Coefficient               0.4739
  Cramer's V                            0.5382


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     39846
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 46152




*/
proc freq data=mme_ucx_flags;
   tables flag_mee_cd19*flag_Constitutive / chisq;
run;

/*
   flag_mee_cd19     flag_Constitutive

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  39847 |    495 |  40342
            |  86.34 |   1.07 |  87.41
            |  98.77 |   1.23 |
            |  91.97 |  17.50 |
   ---------+--------+--------+
          1 |   3477 |   2333 |   5810
            |   7.53 |   5.06 |  12.59
            |  59.85 |  40.15 |
            |   8.03 |  82.50 |
   ---------+--------+--------+
   Total       43324     2828    46152
               93.87     6.13   100.00

              The SAS System         10:39 Monday,

 Statistics for Table of flag_mee_cd19 by flag_Constitutive

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1  13379.4420    <.0001
   Likelihood Ratio Chi-Square    1   8104.5859    <.0001
   Continuity Adj. Chi-Square     1  13372.6752    <.0001
   Mantel-Haenszel Chi-Square     1  13379.1521    <.0001
   Phi Coefficient                       0.5384
   Contingency Coefficient               0.4741
   Cramer's V                            0.5384


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     39847
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 46152

*/
proc freq data=mme_ucx_flags;
   tables flag_mee_tcells*flag_Constitutive / chisq;
run;

/*
 flag_mee_tcells
           flag_Constitutive

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  39968 |    511 |  40479
          |  86.60 |   1.11 |  87.71
          |  98.74 |   1.26 |
          |  92.25 |  18.07 |
 ---------+--------+--------+
        1 |   3356 |   2317 |   5673
          |   7.27 |   5.02 |  12.29
          |  59.16 |  40.84 |
          |   7.75 |  81.93 |
 ---------+--------+--------+
 Total       43324     2828    46152
             93.87     6.13   100.00

Statistics for Table of flag_mee_tcells by flag_Constitutive

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1  13551.3076    <.0001
   Likelihood Ratio Chi-Square    1   8115.7507    <.0001
   Continuity Adj. Chi-Square     1  13544.4275    <.0001
   Mantel-Haenszel Chi-Square     1  13551.0140    <.0001
   Phi Coefficient                       0.5419
   Contingency Coefficient               0.4764
   Cramer's V                            0.5419


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     39968
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 46152


*/
proc freq data=mme_ucx_flags;
   tables flag_mee_all*flag_Constitutive / chisq;
run;


/*
 Table of flag_mee_all by flag_Constitutive

    flag_mee_all     flag_Constitutive

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  40140 |    521 |  40661
             |  86.97 |   1.13 |  88.10
             |  98.72 |   1.28 |
             |  92.65 |  18.42 |
    ---------+--------+--------+
           1 |   3184 |   2307 |   5491
             |   6.90 |   5.00 |  11.90
             |  57.99 |  42.01 |
             |   7.35 |  81.58 |
    ---------+--------+--------+
    Total       43324     2828    46152
                93.87     6.13   100.00


 Statistics for Table of flag_mee_all by flag_Constitutive

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1  13954.1146    <.0001
   Likelihood Ratio Chi-Square    1   8225.6623    <.0001
   Continuity Adj. Chi-Square     1  13947.0341    <.0001
   Mantel-Haenszel Chi-Square     1  13953.8123    <.0001
   Phi Coefficient                       0.5499
   Contingency Coefficient               0.4818
   Cramer's V                            0.5499


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     40140
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 46152


Summary:

While ~50% of MEEs are constitutive, >80% of constitutively present exons are MEEs vs <10% of alternative/unique exons being MEE. This makes sense as if the exon is constitutive, then it is present in all transcripts for a gene and thus should have the highest abundance.

This also makes sense from a functional/evolutionary view: the constitutive exons are those that
are probably more likely to contain functional domains compared to unique/common/alternative exons, in that the functional pieces that define the gene are probably more likely to be conserved AND present in most/all isoforms of a gene. 

*/
