/* Need to determine what percentage of exons were detected

Do for: all genes, autoimmune genes, T1D genes
Then we want the proportion of detected autoimmune among all genes, and proportion of detected T1D among autoimmune */

/* Set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Get detected counts for fusions */

data flag_fusions_on;
   set con.fusions_on_gt_apn0;
   if flag_cd19_on=1 or flag_cd4_on=1 or flag_cd8_on=1 then flag_fusion_on=1;
   else flag_fusion_on=0;
   keep fusion_id flag_fusion_on;
run;


data fusion2gene;
   set fus.unique_info_fusions_si;
   keep fusion_id gene_id flag_multigene;
run;

proc sort data=flag_fusions_on;
   by fusion_id;
proc sort data=fusion2gene;
   by fusion_id;
run;

data flag_fusions_on_w_gene oops;
   merge flag_fusions_on (in=in1) fusion2gene (in=in2);
   by fusion_id;
   if in1 and in2 then output flag_fusions_on_w_gene;
   else output oops;
run;

/* Get immunogene flags */

data immunoflags;
   set con.immunogene_flags;
run;

proc sort data=immunoflags;
   by gene_id;
proc sort data=flag_fusions_on_w_gene;
   by gene_id;
run;

data fusions_on_w_immunogene;
  merge flag_fusions_on_w_gene (in=in1) immunoflags;
  by gene_id;
  if in1;
run;

/* Set missing immunoflags -- these are for multigene fusions! */

data fusions_on_w_immunogene2;
   set fusions_on_w_immunogene;
   if flag_pseudogene=. then flag_pseudogene_multi=0;
   else flag_pseudogene_multi=flag_pseudogene;

   if flag_immuno_gene=. then flag_immuno_gene_multi=0;
   else flag_immuno_gene_multi=flag_immuno_gene;

   if flag_diabetes_gene=. then flag_diabetes_gene_multi=0;
   else flag_diabetes_gene_multi=flag_diabetes_gene;
run;

/* Get percentages of detection for each scenario */

ods listing;
ods html close;


* all genes;
proc freq data=fusions_on_w_immunogene2;
   tables flag_fusion_on;
run;

/*
                                           Cumulative    Cumulative
flag_fusion_on    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0      138168       41.93        138168        41.93
             1      191389       58.07        329557       100.00

58.07% of all exons are detected
*/


* autoimmune genes;
proc freq data=fusions_on_w_immunogene2;
   where flag_immuno_gene=1;
   tables flag_fusion_on;
run;

/*
                                            Cumulative    Cumulative
 flag_fusion_on    Frequency     Percent     Frequency      Percent
 -------------------------------------------------------------------
              0        5005       23.21          5005        23.21
              1       16562       76.79         21567       100.00

76.79% of all autoimmune exons are detected
*/


* T1D genes;
proc freq data=fusions_on_w_immunogene2;
   where flag_diabetes_gene=1;
   tables flag_fusion_on;
run;


/*
                                            Cumulative    Cumulative
 flag_fusion_on    Frequency     Percent     Frequency      Percent
 -------------------------------------------------------------------
              0          62       16.58            62        16.58
              1         312       83.42           374       100.00

83.42% of all T1D exons are detected
*/


/* Now do cross-tabs : what proportion of detected exons are from autoimmune genes ? */

proc freq data=fusions_on_w_immunogene2;
    tables flag_fusion_on*flag_immuno_gene_multi / chisq;
run;


/*
Table of flag_fusion_on by flag_immuno_gene_multi

       flag_fusion_on
                 flag_immuno_gene_multi

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 | 133163 |   5005 | 138168
                |  40.41 |   1.52 |  41.93
                |  96.38 |   3.62 |
                |  43.24 |  23.21 |
       ---------+--------+--------+
              1 | 174827 |  16562 | 191389
                |  53.05 |   5.03 |  58.07
                |  91.35 |   8.65 |
                |  56.76 |  76.79 |
       ---------+--------+--------+
       Total      307990    21567   329557
                   93.46     6.54   100.00


 Statistics for Table of flag_fusion_on by flag_immuno_gene_multi

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1   3320.9935    <.0001
      Likelihood Ratio Chi-Square    1   3550.7816    <.0001
      Continuity Adj. Chi-Square     1   3320.1709    <.0001
      Mantel-Haenszel Chi-Square     1   3320.9834    <.0001
      Phi Coefficient                       0.1004
      Contingency Coefficient               0.0999
      Cramer's V                            0.1004


                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)    133163
                Left-sided Pr <= F          1.0000
                Right-sided Pr >= F         <.0001

                Table Probability (P)       <.0001
                Two-sided Pr <= P           <.0001

                       Sample Size = 329557


Summary:
16562 of 174827 detected fusions are from autoimmune genes (8.65%)
76.79% autoimmune exons detected vs 56.76% non-autoimmune exons

*/

/* what proportion of detected autoimmune exons are from T1D genes ? */

proc freq data=fusions_on_w_immunogene2;
    where flag_immuno_gene=1;
    tables flag_fusion_on*flag_diabetes_gene / chisq;
run;



/*
 Table of flag_fusion_on by flag_diabetes_gene

      flag_fusion_on
                flag_diabetes_gene

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 |   4943 |     62 |   5005
               |  22.92 |   0.29 |  23.21
               |  98.76 |   1.24 |
               |  23.32 |  16.58 |
      ---------+--------+--------+
             1 |  16250 |    312 |  16562
               |  75.35 |   1.45 |  76.79
               |  98.12 |   1.88 |
               |  76.68 |  83.42 |
      ---------+--------+--------+
      Total       21193      374    21567
                  98.27     1.73   100.00

                 The SAS System        14:16 Mond

Statistics for Table of flag_fusion_on by flag_diabetes_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      9.3855    0.0022
   Likelihood Ratio Chi-Square    1     10.1174    0.0015
   Continuity Adj. Chi-Square     1      9.0107    0.0027
   Mantel-Haenszel Chi-Square     1      9.3850    0.0022
   Phi Coefficient                       0.0209
   Contingency Coefficient               0.0209
   Cramer's V                            0.0209


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      4943
             Left-sided Pr <= F          0.9994
             Right-sided Pr >= F         0.0010

             Table Probability (P)       0.0004
             Two-sided Pr <= P           0.0020

                    Sample Size = 21567


Summary:
312 of 16562 detected fusions are from T1D genes (1.88%)
83.42% T1D exons detected vs 76.68% non-T1D autoimmune exons

*/


/* Overall summary:

58.07% of all exons are detected
76.79% of all autoimmune exons are detected
83.42% of all T1D exons are detected

16562 of 174827 detected fusions are from autoimmune genes (8.65%)
76.79% autoimmune exons detected vs 56.76% non-autoimmune exons

312 of 16562 detected fusions are from T1D genes (1.88%)
83.42% T1D exons detected vs 76.68% non-T1D autoimmune exons


*/

