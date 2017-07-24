
ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname event '!MCLAB/event_analysis/sas_data';

/* I want to see if the set of Event-only QDS/DD genes are enriched for autoimmune genes or T1D genes */

/* first, what are the genes in common between MISO and EA? */

data qds_dd_no_miso;
   set event.t1d_miso2qds_comparison2;
   where flag_in_miso=1 and flag_in_qds=1;
   if flag_miso_gene_event_dd=1 then delete;
   if flag_sig_bf5=1 then delete;
run;

data t1d_genes ai_genes;
  set con.immunobase_gene_flags;
  if flag_immuno_gene=1 then output ai_genes;
  if flag_immunobase_diabetes_gene=1 then output t1d_genes;
  keep gene_id;
run;

proc sort data=qds_dd_no_miso;
  by gene_id;
proc sort data=t1d_genes nodup;
  by gene_id;
proc sort data=ai_genes nodup;
  by gene_id;
run;

data miso2qds_immuno;
  merge qds_dd_no_miso (in=in1) t1d_genes (in=in2) ai_genes (in=in3);
  by gene_id;
  if in2 then flag_diabetes_gene=1; else flag_diabetes_gene=0;
  if in3 then flag_immuno_gene=1; else flag_immuno_gene=0;
  if in1 then output;
run;

data flag_as;
   set miso2qds_immuno;
   if flag_cell_by_fus_fdr05=1 or flag_gene_exon_dd=1 then flag_gene_as=1;
   else flag_gene_as=0;
run;

proc freq data=flag_as;
   where flag_miso_testable=1;
   tables flag_gene_as*flag_diabetes_gene flag_gene_as*flag_immuno_gene / chisq;
run;
/*
    flag_gene_as     flag_diabetes_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |    281 |     22 |    303
             |  33.86 |   2.65 |  36.51
             |  92.74 |   7.26 |
             |  35.17 |  70.97 |
    ---------+--------+--------+
           1 |    518 |      9 |    527
             |  62.41 |   1.08 |  63.49
             |  98.29 |   1.71 |
             |  64.83 |  29.03 |
    ---------+--------+--------+
    Total         799       31      830
                96.27     3.73   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1     16.4995    <.0001
 Likelihood Ratio Chi-Square    1     15.7828    <.0001
 Continuity Adj. Chi-Square     1     14.9912    0.0001
 Mantel-Haenszel Chi-Square     1     16.4796    <.0001
 Phi Coefficient                      -0.1410
 Contingency Coefficient               0.1396
 Cramer's V                           -0.1410


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)       281
           Left-sided Pr <= F          <.0001
           Right-sided Pr >= F         1.0000

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                   Sample Size = 830

                     The SAS System             08:46 Monday,

                   The FREQ Procedure

  flag_gene_as     flag_immuno_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |    255 |     48 |    303
           |  30.72 |   5.78 |  36.51
           |  84.16 |  15.84 |
           |  35.71 |  41.38 |
  ---------+--------+--------+
         1 |    459 |     68 |    527
           |  55.30 |   8.19 |  63.49
           |  87.10 |  12.90 |
           |  64.29 |  58.62 |
  ---------+--------+--------+
  Total         714      116      830
              86.02    13.98   100.00


 Statistics for Table of flag_gene_as by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      1.3816    0.2398
  Likelihood Ratio Chi-Square    1      1.3623    0.2431
  Continuity Adj. Chi-Square     1      1.1480    0.2840
  Mantel-Haenszel Chi-Square     1      1.3799    0.2401
  Phi Coefficient                      -0.0408
  Contingency Coefficient               0.0408
  Cramer's V                           -0.0408


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)       255
            Left-sided Pr <= F          0.1422
            Right-sided Pr >= F         0.8990

            Table Probability (P)       0.0411
            Two-sided Pr <= P           0.2535

                    Sample Size = 830

*/


/* Get T1D spliced */

proc print data=flag_as (keep=gene_id flag_miso_testable flag_diabetes_Gene flag_gene_as);
   where  flag_diabetes_gene=1 and flag_gene_as=1;
run;

/*

(of just the testable that are not also DD/DS in MISO)
 C6orf48andSNORD48andSNORD52
 FLOT1
 IKZF1
 KEAP1
 NRM
 PSMD3
 RBM17
 SLC10A3
 TM9SF2


including non testable:
 AIF1
 C6orf48andSNORD48andSNORD52
 CCDC101
 CHST10
 DKC1andSNORA36AandSNORA56
 FLOT1
 IKZF1
 KEAP1
 LIMS3-LOC440895andGPAA1P1
 NRM
 OSM
 PSMD3
 RBM17
 SLC10A3
 TM9SF2


We see IKZF1 and IKZF3, members of the Ikaros family that are alternatively spliced in lymphocytes (refs)
IL7R, which is also AS in lymphocytes, (soluble vs membrane-bound isoforms)
CLEC16A which is thought to be altnernatively spliced in autoimmunity
ERBB3, which expresses different isoforms in different cell lines (check if these are immune cells)
ICAM3, alnteriavely spliced in Chron's disease
PTPN2, AS in autoimmunity

*/

