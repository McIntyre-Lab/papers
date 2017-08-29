ods listing; ods html;

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/concannon/useful_human_data/aceview_hg19/fusions/sas_data';

/* Need to check:

Of genes expressed in 2 or more cell types, are events in general from AI/T1D genes more cell type specific?
	- all events
	- unannotated junctions
	- IR events
	- exon skipping events

Comparisons:
AI vs all other
T1D vs non-AI
T1D vs non-T1D AI
*/

/* Merge cleaned event data with gene lists */

data events_w_flags;
   set splicing.flag_splicing_by_gene_dtct_v2;
run;

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=ai nodup;
   by gene_id;
proc sort data=t1d nodup;
   by gene_id;
proc sort data=events_w_flags;
   by gene_id;
run;

data events_w_flags_immuno;
   merge events_w_flags (in=in1) ai (in=in2) t1d (in=in3);
   by gene_id;
   if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
   if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in1 then output;
run;

/* Flag if event is DD/DE: for this I am classifying any event that is quantitatively different,
   or is detected in one or two cell types as DD/DE */

data flag_event_de_dd;
  set events_w_flags_immuno;
  if sum(flag_cd4_gene_on,flag_cd8_gene_on,flag_cd19_gene_on) > 1 then do;
      if sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) < 3
      and sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) > 0 then flag_event_cell_de_dd=1;
      else if sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) = 3 and flag_anova_fdr_05=1 then flag_event_cell_de_dd=1;
      else if sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) = 3 and flag_anova_fdr_05=0 then flag_event_cell_de_dd=0;
      else flag_event_cell_de_dd=.;
      end;
  else flag_event_cell_de_dd=. ;
run;

/* Compare gene lists */

proc freq data=flag_event_de_dd;
   tables flag_event_cell_de_dd*flag_immuno_gene / chisq ;
run;

/*
 flag_event_cell_de_dd
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  67785 |   8686 |  76471
          |  32.37 |   4.15 |  36.52
          |  88.64 |  11.36 |
          |  36.86 |  34.10 |
 ---------+--------+--------+
        1 | 116126 |  16784 | 132910
          |  55.46 |   8.02 |  63.48
          |  87.37 |  12.63 |
          |  63.14 |  65.90 |
 ---------+--------+--------+
 Total      183911    25470   209381
             87.84    12.16   100.00

Statistics for Table of flag_event_cell_de_dd by flag_immuno_gene

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1     73.2230    <.0001
      Likelihood Ratio Chi-Square    1     73.8513    <.0001
      Continuity Adj. Chi-Square     1     73.1043    <.0001
      Mantel-Haenszel Chi-Square     1     73.2227    <.0001
      Phi Coefficient                       0.0187
      Contingency Coefficient               0.0187
      Cramer's V                            0.0187


                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)     67785
                Left-sided Pr <= F          1.0000
                Right-sided Pr >= F         <.0001

                Table Probability (P)       <.0001
                Two-sided Pr <= P           <.0001

*/

proc freq data=flag_event_de_dd;
   where flag_immuno_gene=1;
   tables flag_event_cell_de_dd*flag_immunobase_diabetes_gene / chisq ;
run;


/*
  flag_event_cell_de_dd
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   6467 |   2219 |   8686
          |  25.39 |   8.71 |  34.10
          |  74.45 |  25.55 |
          |  32.98 |  37.86 |
 ---------+--------+--------+
        1 |  13142 |   3642 |  16784
          |  51.60 |  14.30 |  65.90
          |  78.30 |  21.70 |
          |  67.02 |  62.14 |
 ---------+--------+--------+
 Total       19609     5861    25470
             76.99    23.01   100.00


    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     47.8301    <.0001
    Likelihood Ratio Chi-Square    1     47.2905    <.0001
    Continuity Adj. Chi-Square     1     47.6132    <.0001
    Mantel-Haenszel Chi-Square     1     47.8282    <.0001
    Phi Coefficient                      -0.0433
    Contingency Coefficient               0.0433
    Cramer's V                           -0.0433


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      6467
              Left-sided Pr <= F          <.0001
              Right-sided Pr >= F         1.0000

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001
*/


/* Unannotated junctions */

proc freq data=flag_event_de_dd;
   where flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_event_cell_de_dd*flag_immuno_gene / chisq ;
run;

proc freq data=flag_event_de_dd;
   where flag_immuno_gene=1 and flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_event_cell_de_dd*flag_immunobase_diabetes_gene / chisq ;
run;

/* AUTOIMMUNE GENES vs ALL OTHERS:
Table of flag_event_cell_de_dd by flag_immuno_gene

       flag_event_cell_de_dd
                 flag_immuno_gene

       Frequency|
       Percent  |
       Row Pct  |
       Col Pct  |       0|       1|  Total
       ---------+--------+--------+
              0 |  12277 |   1479 |  13756
                |  31.72 |   3.82 |  35.54
                |  89.25 |  10.75 |
                |  36.02 |  32.05 |
       ---------+--------+--------+
              1 |  21810 |   3135 |  24945
                |  56.36 |   8.10 |  64.46
                |  87.43 |  12.57 |
                |  63.98 |  67.95 |
       ---------+--------+--------+
       Total       34087     4614    38701
                   88.08    11.92   100.00

              Frequency Missing = 16

                  The SAS System              10:49 Frid

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     27.8453    <.0001
   Likelihood Ratio Chi-Square    1     28.2299    <.0001
   Continuity Adj. Chi-Square     1     27.6726    <.0001
   Mantel-Haenszel Chi-Square     1     27.8446    <.0001
   Phi Coefficient                       0.0268
   Contingency Coefficient               0.0268
   Cramer's V                            0.0268


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     12277
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

DIABETES GENES vs NON-T1D AUTOIMMUNE:
 flag_event_cell_de_dd
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   1088 |    391 |   1479
          |  23.58 |   8.47 |  32.05
          |  73.56 |  26.44 |
          |  30.88 |  35.84 |
 ---------+--------+--------+
        1 |   2435 |    700 |   3135
          |  52.77 |  15.17 |  67.95
          |  77.67 |  22.33 |
          |  69.12 |  64.16 |
 ---------+--------+--------+
 Total        3523     1091     4614
             76.35    23.65   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      9.3941    0.0022
 Likelihood Ratio Chi-Square    1      9.2715    0.0023
 Continuity Adj. Chi-Square     1      9.1679    0.0025
 Mantel-Haenszel Chi-Square     1      9.3921    0.0022
 Phi Coefficient                      -0.0451
 Contingency Coefficient               0.0451
 Cramer's V                           -0.0451


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)      1088
           Left-sided Pr <= F          0.0013
           Right-sided Pr >= F         0.9990

           Table Probability (P)       0.0003
           Two-sided Pr <= P           0.0023


*/

/* IR events */

proc freq data=flag_event_de_dd;
   where flag_intron_retention=1;
   tables flag_event_cell_de_dd*flag_immuno_gene / chisq ;
run;

proc freq data=flag_event_de_dd;
   where flag_immuno_gene=1 and flag_intron_retention=1;
   tables flag_event_cell_de_dd*flag_immunobase_diabetes_gene / chisq ;
run;


/* AUTOIMMUNE GENES vs ALL OTHERS:
 Table of flag_event_cell_de_dd by flag_immuno_gene

        flag_event_cell_de_dd
                  flag_immuno_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |  14969 |   1953 |  16922
                 |  36.92 |   4.82 |  41.73
                 |  88.46 |  11.54 |
                 |  41.88 |  40.66 |
        ---------+--------+--------+
               1 |  20775 |   2850 |  23625
                 |  51.24 |   7.03 |  58.27
                 |  87.94 |  12.06 |
                 |  58.12 |  59.34 |
        ---------+--------+--------+
        Total       35744     4803    40547
                    88.15    11.85   100.00

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      2.5758    0.1085
  Likelihood Ratio Chi-Square    1      2.5814    0.1081
  Continuity Adj. Chi-Square     1      2.5260    0.1120
  Mantel-Haenszel Chi-Square     1      2.5757    0.1085
  Phi Coefficient                       0.0080
  Contingency Coefficient               0.0080
  Cramer's V                            0.0080


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     14969
            Left-sided Pr <= F          0.9476
            Right-sided Pr >= F         0.0559

            Table Probability (P)       0.0034
            Two-sided Pr <= P           0.1120

DIABETES GENES vs NON-T1D AUTOIMMUNE:
  flag_event_cell_de_dd
            flag_immunobase_diabetes_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   1459 |    494 |   1953
           |  30.38 |  10.29 |  40.66
           |  74.71 |  25.29 |
           |  39.80 |  43.45 |
  ---------+--------+--------+
         1 |   2207 |    643 |   2850
           |  45.95 |  13.39 |  59.34
           |  77.44 |  22.56 |
           |  60.20 |  56.55 |
  ---------+--------+--------+
  Total        3666     1137     4803
              76.33    23.67   100.00

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1      4.7906    0.0286
    Likelihood Ratio Chi-Square    1      4.7691    0.0290
    Continuity Adj. Chi-Square     1      4.6405    0.0312
    Mantel-Haenszel Chi-Square     1      4.7896    0.0286
    Phi Coefficient                      -0.0316
    Contingency Coefficient               0.0316
    Cramer's V                           -0.0316


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      1459
              Left-sided Pr <= F          0.0158
              Right-sided Pr >= F         0.9868

              Table Probability (P)       0.0025
              Two-sided Pr <= P           0.0295

                 Effective Sample Size = 4803

*/

/* Exon skipping events */

proc freq data=flag_event_de_dd;
   where flag_exonskip=1;
   tables flag_event_cell_de_dd*flag_immuno_gene / chisq ;
run;


proc freq data=flag_event_de_dd;
   where flag_immuno_gene=1 and flag_exonskip=1;
   tables flag_event_cell_de_dd*flag_immunobase_diabetes_gene / chisq ;
run;


/* AUTOIMMUNE GENES vs ALL OTHERS:

 flag_event_cell_de_dd
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  14768 |   1639 |  16407
          |  32.06 |   3.56 |  35.62
          |  90.01 |   9.99 |
          |  36.08 |  31.94 |
 ---------+--------+--------+
        1 |  26162 |   3492 |  29654
          |  56.80 |   7.58 |  64.38
          |  88.22 |  11.78 |
          |  63.92 |  68.06 |
 ---------+--------+--------+
 Total       40930     5131    46061
             88.86    11.14   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1     34.0448    <.0001
 Likelihood Ratio Chi-Square    1     34.5451    <.0001
 Continuity Adj. Chi-Square     1     33.8646    <.0001
 Mantel-Haenszel Chi-Square     1     34.0440    <.0001
 Phi Coefficient                       0.0272
 Contingency Coefficient               0.0272
 Cramer's V                            0.0272


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     14768
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

             Effective Sample Size = 46061
                 Frequency Missing = 20


DIABETES GENES vs NON-T1D AUTOIMMUNE:
    flag_event_cell_de_dd
              flag_immunobase_diabetes_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |   1260 |    379 |   1639
             |  24.56 |   7.39 |  31.94
             |  76.88 |  23.12 |
             |  31.25 |  34.49 |
    ---------+--------+--------+
           1 |   2772 |    720 |   3492
             |  54.02 |  14.03 |  68.06
             |  79.38 |  20.62 |
             |  68.75 |  65.51 |
    ---------+--------+--------+
    Total        4032     1099     5131
                78.58    21.42   100.00

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      4.1596    0.0414
   Likelihood Ratio Chi-Square    1      4.1196    0.0424
   Continuity Adj. Chi-Square     1      4.0121    0.0452
   Mantel-Haenszel Chi-Square     1      4.1588    0.0414
   Phi Coefficient                      -0.0285
   Contingency Coefficient               0.0285
   Cramer's V                           -0.0285


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)      1260
             Left-sided Pr <= F          0.0230
             Right-sided Pr >= F         0.9807

             Table Probability (P)       0.0037
             Two-sided Pr <= P           0.0447

*/
