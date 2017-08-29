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

/* Flag if event is cell type specific */

data flag_event_spec;
  set events_w_flags_immuno;
  if sum(flag_cd4_gene_on,flag_cd8_gene_on,flag_cd19_gene_on) > 1 then do;
      if sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) = 1 then flag_event_cell_specific=1;
      else if sum(flag_cd4_on,flag_cd8_on,flag_cd19_on) > 1 then flag_event_cell_specific=0;
      end;
  else flag_event_cell_specific=. ;
  /* Set up my non-AI flag so I can compare between non-AI genes and T1D genes */
  if flag_immuno_gene=0 and flag_immunobase_diabetes_gene=0 then flag_diabetes_vs_others=0;
  else if flag_immuno_gene=1 and flag_immunobase_diabetes_gene=1 then flag_diabetes_vs_others=1;
  else flag_diabetes_vs_others=.;
run;

/* Compare gene lists */

proc freq data=flag_event_spec;
   tables flag_event_cell_specific*flag_immuno_gene / chisq ;
run;

/*
flag_event_cell_specific
          flag_immuno_gene

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 | 165721 |  22988 | 188709
         |  79.15 |  10.98 |  90.13
         |  87.82 |  12.18 |
         |  90.11 |  90.26 |
---------+--------+--------+
       1 |  18190 |   2482 |  20672
         |   8.69 |   1.19 |   9.87
         |  87.99 |  12.01 |
         |   9.89 |   9.74 |
---------+--------+--------+
Total      183911    25470   209381
            87.84    12.16   100.00

      Frequency Missing = 506

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1      0.5349    0.4646
   Likelihood Ratio Chi-Square    1      0.5367    0.4638
   Continuity Adj. Chi-Square     1      0.5186    0.4714
   Mantel-Haenszel Chi-Square     1      0.5349    0.4646
   Phi Coefficient                      -0.0016
   Contingency Coefficient               0.0016
   Cramer's V                           -0.0016


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    165721
             Left-sided Pr <= F          0.2361
             Right-sided Pr >= F         0.7708

             Table Probability (P)       0.0069
             Two-sided Pr <= P           0.4732


*/

proc freq data=flag_event_spec;
   where flag_immuno_gene=1;
   tables flag_event_cell_specific*flag_immunobase_diabetes_gene / chisq ;
run;


/*
 flag_event_cell_specific
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  17705 |   5283 |  22988
          |  69.51 |  20.74 |  90.26
          |  77.02 |  22.98 |
          |  90.29 |  90.14 |
 ---------+--------+--------+
        1 |   1904 |    578 |   2482
          |   7.48 |   2.27 |   9.74
          |  76.71 |  23.29 |
          |   9.71 |   9.86 |
 ---------+--------+--------+
 Total       19609     5861    25470
             76.99    23.01   100.00

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1      0.1185    0.7307
     Likelihood Ratio Chi-Square    1      0.1182    0.7310
     Continuity Adj. Chi-Square     1      0.1018    0.7496
     Mantel-Haenszel Chi-Square     1      0.1185    0.7307
     Phi Coefficient                       0.0022
     Contingency Coefficient               0.0022
     Cramer's V                            0.0022


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)     17705
               Left-sided Pr <= F          0.6452
               Right-sided Pr >= F         0.3736

               Table Probability (P)       0.0188
               Two-sided Pr <= P           0.7254

*/


/* Unannotated junctions */

proc freq data=flag_event_spec;
   where flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_event_cell_specific*flag_immuno_gene / chisq ;
run;

proc freq data=flag_event_spec;
   where flag_immuno_gene=1 and flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_event_cell_specific*flag_immunobase_diabetes_gene / chisq ;
run;

/* AUTOIMMUNE GENES vs ALL OTHERS:
 flag_event_cell_specific
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  30013 |   4068 |  34081
          |  77.55 |  10.51 |  88.06
          |  88.06 |  11.94 |
          |  88.05 |  88.17 |
 ---------+--------+--------+
        1 |   4074 |    546 |   4620
          |  10.53 |   1.41 |  11.94
          |  88.18 |  11.82 |
          |  11.95 |  11.83 |
 ---------+--------+--------+
 Total       34087     4614    38701
             88.08    11.92   100.00

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      0.0540    0.8162
  Likelihood Ratio Chi-Square    1      0.0541    0.8160
  Continuity Adj. Chi-Square     1      0.0434    0.8350
  Mantel-Haenszel Chi-Square     1      0.0540    0.8162
  Phi Coefficient                      -0.0012
  Contingency Coefficient               0.0012
  Cramer's V                           -0.0012


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     30013
            Left-sided Pr <= F          0.4193
            Right-sided Pr >= F         0.5996

            Table Probability (P)       0.0188
            Two-sided Pr <= P           0.8277


DIABETES GENES vs NON-T1D AUTOIMMUNE:


   flag_event_cell_specific
             flag_immunobase_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |   3118 |    950 |   4068
            |  67.58 |  20.59 |  88.17
            |  76.65 |  23.35 |
            |  88.50 |  87.08 |
   ---------+--------+--------+
          1 |    405 |    141 |    546
            |   8.78 |   3.06 |  11.83
            |  74.18 |  25.82 |
            |  11.50 |  12.92 |
   ---------+--------+--------+
   Total        3523     1091     4614
               76.35    23.65   100.00


  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      1.6283    0.2019
  Likelihood Ratio Chi-Square    1      1.5999    0.2059
  Continuity Adj. Chi-Square     1      1.4943    0.2216
  Mantel-Haenszel Chi-Square     1      1.6279    0.2020
  Phi Coefficient                       0.0188
  Contingency Coefficient               0.0188
  Cramer's V                            0.0188


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)      3118
            Left-sided Pr <= F          0.9073
            Right-sided Pr >= F         0.1114

            Table Probability (P)       0.0187
            Two-sided Pr <= P           0.2172

*/

/* IR events */

proc freq data=flag_event_spec;
   where flag_intron_retention=1;
   tables flag_event_cell_specific*flag_immuno_gene / chisq ;
run;

proc freq data=flag_event_spec;
   where flag_immuno_gene=1 and flag_intron_retention=1;
   tables flag_event_cell_specific*flag_immunobase_diabetes_gene / chisq ;
run;


/* AUTOIMMUNE GENES vs ALL OTHERS:

    flag_event_cell_specific
              flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  30011 |   3945 |  33956
             |  74.02 |   9.73 |  83.74
             |  88.38 |  11.62 |
             |  83.96 |  82.14 |
    ---------+--------+--------+
           1 |   5733 |    858 |   6591
             |  14.14 |   2.12 |  16.26
             |  86.98 |  13.02 |
             |  16.04 |  17.86 |
    ---------+--------+--------+
    Total       35744     4803    40547
                88.15    11.85   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1     10.3568    0.0013
 Likelihood Ratio Chi-Square    1     10.1304    0.0015
 Continuity Adj. Chi-Square     1     10.2232    0.0014
 Mantel-Haenszel Chi-Square     1     10.3566    0.0013
 Phi Coefficient                       0.0160
 Contingency Coefficient               0.0160
 Cramer's V                            0.0160


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     30011
           Left-sided Pr <= F          0.9993
           Right-sided Pr >= F         0.0008

           Table Probability (P)       0.0001
           Two-sided Pr <= P           0.0014

DIABETES GENES vs NON-T1D AUTOIMMUNE:
   flag_event_cell_specific
             flag_immunobase_diabetes_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |   3006 |    939 |   3945
            |  62.59 |  19.55 |  82.14
            |  76.20 |  23.80 |
            |  82.00 |  82.59 |
   ---------+--------+--------+
          1 |    660 |    198 |    858
            |  13.74 |   4.12 |  17.86
            |  76.92 |  23.08 |
            |  18.00 |  17.41 |
   ---------+--------+--------+
   Total        3666     1137     4803
               76.33    23.67   100.00

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      0.2052    0.6505
  Likelihood Ratio Chi-Square    1      0.2062    0.6498
  Continuity Adj. Chi-Square     1      0.1670    0.6828
  Mantel-Haenszel Chi-Square     1      0.2052    0.6506
  Phi Coefficient                      -0.0065
  Contingency Coefficient               0.0065
  Cramer's V                           -0.0065


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)      3006
            Left-sided Pr <= F          0.3429
            Right-sided Pr >= F         0.6892

            Table Probability (P)       0.0321
            Two-sided Pr <= P           0.6901

*/

/* Exon skipping events */

proc freq data=flag_event_spec;
   where flag_exonskip=1;
   tables flag_event_cell_specific*flag_immuno_gene / chisq ;
run;


proc freq data=flag_event_spec;
   where flag_immuno_gene=1 and flag_exonskip=1;
   tables flag_event_cell_specific*flag_immunobase_diabetes_gene / chisq ;
run;


/* AUTOIMMUNE GENES vs ALL OTHERS:
 flag_event_cell_specific
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  36153 |   4493 |  40646
          |  78.49 |   9.75 |  88.24
          |  88.95 |  11.05 |
          |  88.33 |  87.57 |
 ---------+--------+--------+
        1 |   4777 |    638 |   5415
          |  10.37 |   1.39 |  11.76
          |  88.22 |  11.78 |
          |  11.67 |  12.43 |
 ---------+--------+--------+
 Total       40930     5131    46061
             88.86    11.14   100.00

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      2.5592    0.1097
  Likelihood Ratio Chi-Square    1      2.5229    0.1122
  Continuity Adj. Chi-Square     1      2.4861    0.1149
  Mantel-Haenszel Chi-Square     1      2.5591    0.1097
  Phi Coefficient                       0.0075
  Contingency Coefficient               0.0075
  Cramer's V                            0.0075


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     36153
            Left-sided Pr <= F          0.9469
            Right-sided Pr >= F         0.0582

            Table Probability (P)       0.0051
            Two-sided Pr <= P           0.1126


DIABETES GENES vs NON-T1D AUTOIMMUNE:
  flag_event_cell_specific
            flag_immunobase_diabetes_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   3548 |    945 |   4493
           |  69.15 |  18.42 |  87.57
           |  78.97 |  21.03 |
           |  88.00 |  85.99 |
  ---------+--------+--------+
         1 |    484 |    154 |    638
           |   9.43 |   3.00 |  12.43
           |  75.86 |  24.14 |
           |  12.00 |  14.01 |
  ---------+--------+--------+
  Total        4032     1099     5131
              78.58    21.42   100.00

      Statistic                     DF       Value      Prob
      ------------------------------------------------------
      Chi-Square                     1      3.2005    0.0736
      Likelihood Ratio Chi-Square    1      3.1217    0.0773
      Continuity Adj. Chi-Square     1      3.0187    0.0823
      Mantel-Haenszel Chi-Square     1      3.1999    0.0736
      Phi Coefficient                       0.0250
      Contingency Coefficient               0.0250
      Cramer's V                            0.0250


                       Fisher's Exact Test
                ----------------------------------
                Cell (1,1) Frequency (F)      3548
                Left-sided Pr <= F          0.9659
                Right-sided Pr >= F         0.0424

                Table Probability (P)       0.0083
                Two-sided Pr <= P           0.0795

*/
