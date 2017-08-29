/* Set libraries */

ods listing; ods html close;

libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname fus '!MCLAB/useful_human_data/aceview_hg19/fusions/sas_data';
libname av '!MCLAB/useful_human_data/aceview_hg19/sas_data';

/* Re-import immunogene list and T1D list */

proc import datafile='!MCLAB/useful_human_data/aceview_hg19/downloaded_files/immunogenes_hg19.csv'
     out=immunobase_genes dbms=csv replace; guessingrows=15791;
run; quit;

proc import datafile='!MCLAB/useful_human_data/aceview_hg19/downloaded_files/immunobase_hg19_t1d_genes.csv'
     out=diabetes_genes dbms=csv replace; guessingrows=1026;
run; quit;


/* Get list of aceview genes with entrezID */

data diabetes_genes2;
   set diabetes_genes;
   if EntrezGene_ID_s_=. then delete; *drop genes without an entrezID, as there is no clean way to merge them;
   rename EntrezGene_ID_s_=entrez_id hgnc_symbol=hgnc_id;
   drop biotype cand_gene in_region;
run;

data immunobase_genes2;
   set immunobase_genes;
   if entrez_id="NA" then delete;  *drop genes without an entrezID, as there is no clean way to merge them;
   entrez_id2=entrez_id * 1; *convert entrez_id from string to numeric;
   drop gene_type entrez_id;
   rename entrez_id2=entrez_id;
run;


proc sort data=av.aceview2entrez;
   by entrez_id;
proc sort data=immunobase_genes2 nodup;
   by entrez_id;
proc sort data=diabetes_genes2  nodup;
   by entrez_id;
run;

data aceview_w_ai no_av;
   merge av.aceview2entrez (in=in1) immunobase_genes2 (in=in2);
   by entrez_id;
   if entrez_id ne . then do;
      if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
      if in1 then output aceview_w_ai; else output no_av; *294 genes not in Aceview;
      end;
   else do;
     flag_immuno_gene=0;
     output aceview_w_ai;
     end;
run;

proc sort data=aceview_w_ai nodup;
   by entrez_id;
run;


data aceview_w_flags no_av;
   merge aceview_w_ai (in=in1) diabetes_genes2 (in=in2);
   by entrez_id;
   if entrez_id ne . then do;
      if in2 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
      if in1 then output aceview_w_flags; else output no_av; *29 genes not in Aceview;
      end;
   else do;
     flag_immunobase_diabetes_gene=0;
     output aceview_w_flags;
     end;
run;
* Some entrez_ids are duplicated so we can sort uniquely to get the final set;

data aceview_w_flags2;
   set aceview_w_flags;
   keep gene_id flag_immuno_gene flag_immunobase_diabetes_gene;
run;

proc sort data=aceview_w_flags2 nodup;
   by gene_id;
run;

*383485 genes in total!;

/* Make permenant */

data con.immunobase_gene_flags;
  set aceview_w_flags2;
run;


proc freq data=con.immunobase_gene_flags;
   tables flag_immuno_gene*flag_immunobase_diabetes_gene;
run;

*1852 autoimmune genes;
*432 diabetes genes;

/* Assign immunogene flags to fusions (makes it easier later!) */

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   keep fusion_id gene_id;
run;

proc sort data=ai nodup;
   by gene_id;
proc sort data=t1d  nodup;
   by gene_id;
proc sort data=fus2gene;
   by gene_id;
run;
   
   
data ai_fus;
   merge fus2gene (in=in1) ai (in=in2);
   by gene_id;
   if in1 and in2;
   keep fusion_id;
run;

data t1d_fus;
   merge fus2gene (in=in1) t1d (in=in2);
   by gene_id;
   if in1 and in2;
   keep fusion_id;
run;

data fusions;
   set fus.unique_info_fusions_si;
   keep fusion_id;
   run;

proc sort data=ai_fus nodup;
  by fusion_id;
proc sort data=t1d_fus nodup;
  by fusion_id;
proc sort data=fusions nodup;
  by fusion_id;
run;

data fusions_w_gene_flags;
  merge fusions (in=in1) ai_fus (in=in2) t1d_fus (in=in3);
  by fusion_id;
  if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
  if in3 then flag_ibase_gene=1; else flag_ibase_gene=0;
  run;
  
  

/* Redo fusion and splicing counts on ALL immunobase T1D genes */
/* use flag_immunobase_diabetes_gene in con.diabetes_gene_list_full */


* fusion detection;
data flag_fusions_on;
   set con.fusions_on_gt_apn0;
   if flag_cd19_on=1 or flag_cd4_on=1 or flag_cd8_on=1 then flag_fusion_on=1;
   else flag_fusion_on=0;
   keep fusion_id flag_fusion_on flag_cd19_on flag_cd4_on flag_cd8_on;
run;


proc sort data=flag_fusions_on;
  by fusion_id;
proc sort data=fusions_w_gene_flags;
  by fusion_id;
run;


data fusions_on_w_flags;
   merge flag_fusions_on (in=in1) fusions_w_gene_flags (in=in2);
   by fusion_id;
   if in1 and in2;
   run;
   
proc freq data=fusions_on_w_flags;
   where flag_immuno_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=ai_fus_counts;
   tables flag_fusion_on;
run;

*18349 of 23539 fusions detected;

proc print data=ai_fus_counts;
run;

/*
  flag_     flag_     flag_
 CD4_on    CD8_on    CD19_on    COUNT

    0         0         0        5190
    1         0         0         176
    0         1         0         234
    0         0         1         494
    1         0         1         110
    1         1         0         815
    0         1         1         165
    1         1         1       16355
*/

proc freq data=fusions_on_w_flags;
   where flag_ibase_gene=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on  / out=t1d_fus_counts;
   tables flag_fusion_on;
run;
   
*3969 of 5002 fusions on;

proc print data=t1d_fus_counts;
run;

/*
 flag_     flag_     flag_
CD4_on    CD8_on    CD19_on    COUNT

   0         0         0        1033
   1         0         0          31
   0         1         0          39
   0         0         1         101
   1         0         1          23
   1         1         0         247
   0         1         1          35
   1         1         1        3493


*/

* fusion DE;


data fus_results_for_enrich;
   set con.results_by_fusion_w_flags;

  /* Scenario 1 - only DE fusion expressed in all three cell types */

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then flag_fusion_de=.;
   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fusion_de=1;
   else flag_fusion_de=0;

  /* Scenario 2 - only DE fusion expressed in all three cell types PLUS fusions expressed in only one cell type */

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then do;
        if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_fusion_dd_1=1;
        else if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_fusion_dd_1=1;
        else if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_fusion_dd_1=1;
        else flag_fusion_dd_1=.;
        end;

   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fusion_dd_1=1;
   else flag_fusion_dd_1=0;

  /* Scenario 3 - only DE fusion expressed in all three cell types PLUS fusions expressed in only one or two cell types */

   if flag_cd4cd8_fdr05=. and flag_cd4cd19_fdr05=. and flag_cd8cd19_fdr05=. then do;
     if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_fusion_dd_2=1;
     else flag_fusion_dd_2=.;
     end;
   else if flag_cd4cd8_fdr05=1 or flag_cd4cd19_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_fusion_dd_2=1;
   else flag_fusion_dd_2=0;
   

   if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_any_on=1;
   else flag_any_on=0;
   drop flag_immuno_gene flag_diabetes_gene;
   
run;

proc sort data=fusions_w_gene_flags;
  by fusion_id;
proc sort data=fus_results_for_enrich;
  by fusion_id;
run;

data fus_results_de;
    merge fus_results_for_enrich (in=in1) fusions_w_gene_flags (in=in2);
    by fusion_id;
    if in1 and in2;
run;


/* Autoimmune genes 
DE=exp all 3 tissues and quantitatively different
DD1 = DE + only detected in one tissue
DD2 = DE + only detected in one or two tissues
*/

proc freq data=fus_results_de;
   where flag_immuno_gene=1 and flag_any_on=1;
   tables flag_fusion_de flag_fusion_dd_1 flag_fusion_dd_2;
   tables flag_cd4cd8_fdr05*flag_cd4cd19_fdr05*flag_cd8cd19_fdr05 / out=ai_de_count;
run;

/* DE = 13262; DD1 = 13438; DD2 = 15256 */

proc print data=ai_de_count;
run;

/*
 flag_       flag_       flag_
cd4cd8_    cd4cd19_    cd8cd19_
 fdr05       fdr05       fdr05     COUNT    PERCENT
.           .           .        1994      .
0           0           0        3093    18.9116
0           0           1         740     4.5246
0           1           0         839     5.1299
0           1           1        6972    42.6292
1           0           0         147     0.8988
1           0           1         647     3.9560
1           1           0         603     3.6869
1           1           1        3314    20.2629
*/

/* Autoimmune genes 
DE=exp all 3 tissues and quantitatively different
DD1 = DE + only detected in one tissue
DD2 = DE + only detected in one or two tissues
*/

proc freq data=fus_results_de;
   where flag_ibase_gene=1 and flag_any_on=1;
   tables flag_fusion_de flag_fusion_dd_1 flag_fusion_dd_2;
   tables flag_cd4cd8_fdr05*flag_cd4cd19_fdr05*flag_cd8cd19_fdr05 / out=t1d_de_count;
run;

/* DE = 2717; DD1 = 2748; DD2 = 3193 */

proc print data=t1d_de_count;
run;

/*
 flag_       flag_       flag_
cd4cd8_    cd4cd19_    cd8cd19_
 fdr05       fdr05       fdr05     COUNT    PERCENT

   .           .           .         476      .
   0           0           0         776    22.2159
   0           0           1         152     4.3516
   0           1           0         187     5.3536
   0           1           1        1347    38.5628
   1           0           0          30     0.8589
   1           0           1         130     3.7217
   1           1           0         146     4.1798
   1           1           1         725    20.7558

*/


/* Redo fusion enrichments:
   T1D vs all autoimmune
   T1D vs all genes */

* Detection enrichment ;

* Autoimmune vs all;

proc freq data=fus_results_de;
   tables flag_any_on*flag_immuno_gene / chisq;
run;

proc freq data=fusions_on_w_flags;
   tables flag_fusion_on*flag_immuno_gene / chisq;
run;


/*  flag_any_on     flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 132978 |   5190 | 138168
          |  40.35 |   1.57 |  41.93
          |  96.24 |   3.76 |
          |  43.45 |  22.05 |
 ---------+--------+--------+
        1 | 173040 |  18349 | 191389
          |  52.51 |   5.57 |  58.07
          |  90.41 |   9.59 |
          |  56.55 |  77.95 |
 ---------+--------+--------+
 Total      306018    23539   329557
             92.86     7.14   100.00

 Statistics for Table of flag_any_on by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1   4113.4257    <.0001
  Likelihood Ratio Chi-Square    1   4424.5417    <.0001
  Continuity Adj. Chi-Square     1   4112.5466    <.0001
  Mantel-Haenszel Chi-Square     1   4113.4132    <.0001
  Phi Coefficient                       0.1117
  Contingency Coefficient               0.1110
  Cramer's V                            0.1117


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)    132978
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 329557

*/

* T1D vs Autoimmune;

proc freq data=fus_results_de;
   where flag_immuno_gene=1;
   tables flag_any_on*flag_ibase_gene / chisq;
run;

proc freq data=fusions_on_w_flags;
   where flag_immuno_gene=1;
   tables flag_fusion_on*flag_ibase_gene / chisq;
run;


/*
 flag_fusion_on     flag_ibase_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   4157 |   1033 |   5190
          |  17.66 |   4.39 |  22.05
          |  80.10 |  19.90 |
          |  22.43 |  20.65 |
 ---------+--------+--------+
        1 |  14380 |   3969 |  18349
          |  61.09 |  16.86 |  77.95
          |  78.37 |  21.63 |
          |  77.57 |  79.35 |
 ---------+--------+--------+
 Total       18537     5002    23539
             78.75    21.25   100.00

Statistics for Table of flag_fusion_on by flag_ibase_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      7.2101    0.0072
  Likelihood Ratio Chi-Square    1      7.2930    0.0069
  Continuity Adj. Chi-Square     1      7.1073    0.0077
  Mantel-Haenszel Chi-Square     1      7.2098    0.0073
  Phi Coefficient                       0.0175
  Contingency Coefficient               0.0175
  Cramer's V                            0.0175


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)      4157
            Left-sided Pr <= F          0.9967
            Right-sided Pr >= F         0.0037

            Table Probability (P)       0.0004
            Two-sided Pr <= P           0.0071

                   Sample Size = 23539
*/

* DE/DD enrichment ;
/* DE enrichment - AI genes */
proc freq data=fus_results_de;
   tables flag_fusion_de*flag_immuno_gene / chisq;
run;

/*
 Table of flag_fusion_de by flag_immuno_gene

     flag_fusion_de     flag_immuno_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |  35461 |   3093 |  38554
              |  21.66 |   1.89 |  23.55
              |  91.98 |   8.02 |
              |  24.06 |  18.91 |
     ---------+--------+--------+
            1 | 111897 |  13262 | 125159
              |  68.35 |   8.10 |  76.45
              |  89.40 |  10.60 |
              |  75.94 |  81.09 |
     ---------+--------+--------+
     Total      147358    16355   163713
                 90.01     9.99   100.00

 Statistics for Table of flag_fusion_de by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    217.1077    <.0001
   Likelihood Ratio Chi-Square    1    227.0016    <.0001
   Continuity Adj. Chi-Square     1    216.8216    <.0001
   Mantel-Haenszel Chi-Square     1    217.1063    <.0001
   Phi Coefficient                       0.0364
   Contingency Coefficient               0.0364
   Cramer's V                            0.0364

        Fisher's Exact Test
 ----------------------------------
 Cell (1,1) Frequency (F)     35461
 Left-sided Pr <= F          1.0000
 Right-sided Pr >= F         <.0001

 Table Probability (P)       <.0001
 Two-sided Pr <= P           <.0001

   Effective Sample Size = 163713

*/

/* DE enrichment - IBase vs AI*/
proc freq data=fus_results_de;
   where flag_fusion_all_on0=1 and flag_immuno_gene=1;
   tables flag_fusion_de*flag_ibase_gene / chisq;
run;

/*
 Table of flag_fusion_de by flag_ibase_gene

    flag_fusion_de     flag_ibase_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |   2317 |    776 |   3093
             |  14.17 |   4.74 |  18.91
             |  74.91 |  25.09 |
             |  18.01 |  22.22 |
    ---------+--------+--------+
           1 |  10545 |   2717 |  13262
             |  64.48 |  16.61 |  81.09
             |  79.51 |  20.49 |
             |  81.99 |  77.78 |
    ---------+--------+--------+
    Total       12862     3493    16355
                78.64    21.36   100.00

Statistics for Table of flag_fusion_de by flag_ibase_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     31.6220    <.0001
  Likelihood Ratio Chi-Square    1     30.6898    <.0001
  Continuity Adj. Chi-Square     1     31.3486    <.0001
  Mantel-Haenszel Chi-Square     1     31.6200    <.0001
  Phi Coefficient                      -0.0440
  Contingency Coefficient               0.0439
  Cramer's V                           -0.0440


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)      2317
            Left-sided Pr <= F          <.0001
            Right-sided Pr >= F         1.0000

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 16355
*/


/* DD2 enrichment - AI genes */
proc freq data=fus_results_de;
   tables flag_fusion_dd_2*flag_immuno_gene / chisq;
run;

/*
 flag_fusion_dd_2
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  35461 |   3093 |  38554
          |  18.53 |   1.62 |  20.14
          |  91.98 |   8.02 |
          |  20.49 |  16.86 |
 ---------+--------+--------+
        1 | 137579 |  15256 | 152835
          |  71.88 |   7.97 |  79.86
          |  90.02 |   9.98 |
          |  79.51 |  83.14 |
 ---------+--------+--------+
 Total      173040    18349   191389
             90.41     9.59   100.00

Statistics for Table of flag_fusion_dd_2 by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1    136.3759    <.0001
   Likelihood Ratio Chi-Square    1    141.8281    <.0001
   Continuity Adj. Chi-Square     1    136.1500    <.0001
   Mantel-Haenszel Chi-Square     1    136.3752    <.0001
   Phi Coefficient                       0.0267
   Contingency Coefficient               0.0267
   Cramer's V                            0.0267


         Fisher's Exact Test
  ----------------------------------
  Cell (1,1) Frequency (F)     35461
  Left-sided Pr <= F          1.0000
  Right-sided Pr >= F         <.0001

  Table Probability (P)       <.0001
  Two-sided Pr <= P           <.0001

    Effective Sample Size = 191389

*/

/* DD2 enrichment - IBase vs AI*/
proc freq data=fus_results_de;
   where flag_fusion_all_on0=1 and flag_immuno_gene=1;
   tables flag_fusion_dd_2*flag_ibase_gene / chisq;
run;

/*
Table of flag_fusion_dd_2 by flag_ibase_gen

    flag_fusion_dd_2
              flag_ibase_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |   2317 |    776 |   3093
             |  14.17 |   4.74 |  18.91
             |  74.91 |  25.09 |
             |  18.01 |  22.22 |
    ---------+--------+--------+
           1 |  10545 |   2717 |  13262
             |  64.48 |  16.61 |  81.09
             |  79.51 |  20.49 |
             |  81.99 |  77.78 |
    ---------+--------+--------+
    Total       12862     3493    16355
                78.64    21.36   100.00

 Statistics for Table of flag_fusion_dd_2 by flag_ibase_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     31.6220    <.0001
    Likelihood Ratio Chi-Square    1     30.6898    <.0001
    Continuity Adj. Chi-Square     1     31.3486    <.0001
    Mantel-Haenszel Chi-Square     1     31.6200    <.0001
    Phi Coefficient                      -0.0440
    Contingency Coefficient               0.0439
    Cramer's V                           -0.0440


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)      2317
              Left-sided Pr <= F          <.0001
              Right-sided Pr >= F         1.0000

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 16355
*/



/* Splicing counts */

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=ai nodup;
   by gene_id;
proc sort data=t1d  nodup;
   by gene_id;
proc sort data=splicing.splicing_results_clean;
   by gene_id;
run;

data splicing_results_clean_w_flags;
  merge splicing.splicing_results_clean (in=in1) t1d (in=in2) ai (in=in3);
   by gene_id;
   if in2 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in3 then flag_immuno_gene=1; else flag_immuno_gene=0;
   if in1 then output;
run;


/* Get counts -- made a macro as the code is virtually the same for both autoimmune and T1D genes*/

%macro counts(flag);

* total splicing events;
proc freq data=splicing_results_clean_w_flags noprint ;
   where &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_total;
run;

* exon skipping events;
proc freq data=splicing_results_clean_w_flags noprint ;
   where flag_exonskip=1 and &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_es;
run;

* alt donor events;
proc freq data=splicing_results_clean_w_flags noprint ;
   where flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_ad;
run;

* alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint ;
   where flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_aa;
run;

* alt donor and alt acceptor events;
proc freq data=splicing_results_clean_w_flags noprint ;
   where flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_ada;
run;

* intron retention events;
proc freq data=splicing_results_clean_w_flags noprint ;
   where flag_intron_retention=1 and &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_ir;
run;

* annotated junctions;
proc freq data=splicing_results_clean_w_flags noprint ;
   where flag_junction_annotated=1 and &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_annot;
run;

* unannotated junctions;
proc freq data=splicing_results_clean_w_flags noprint;
   where flag_junction_annotated=0 and flag_intron_retention=0 and &flag.=1;
   tables flag_cd4_on*flag_cd8_on*flag_cd19_on / out=splicing_unannot;
run;


/* format datasets and merge so I can print out one table of counts */

%macro format(dataset,class);

data &dataset._2;
   length cell $8.;
   set &dataset.;
   if flag_cd4_on=. and flag_cd8_on=. and flag_cd19_on=. then cell="missing";
   if flag_cd4_on=0 and flag_cd8_on=0 and flag_cd19_on=0 then cell="off";
   if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then cell="cd4";
   if flag_cd4_on=0 and flag_cd8_on=1 and flag_cd19_on=0 then cell="cd8";
   if flag_cd4_on=0 and flag_cd8_on=0 and flag_cd19_on=1 then cell="cd19";
   if flag_cd4_on=1 and flag_cd8_on=1 and flag_cd19_on=0 then cell="cd4_8";
   if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=1 then cell="cd4_19";
   if flag_cd4_on=0 and flag_cd8_on=1 and flag_cd19_on=1 then cell="cd8_19";
   if flag_cd4_on=1 and flag_cd8_on=1 and flag_cd19_on=1 then cell="cd4_8_19";
   keep cell count;
   rename count=&class._count;
run;

proc sort data=&dataset._2;
  by cell;
run;

%mend;

%format(splicing_total,total);
%format(splicing_es,ES);
%format(splicing_ad,AD);
%format(splicing_aa,AA);
%format(splicing_ada,ADA);
%format(splicing_ir,IR);
%format(splicing_annot,ANNOT);
%format(splicing_unannot,UNANNOT);

data all_counts;
  merge splicing_total_2 splicing_annot_2 splicing_unannot_2 splicing_ir_2
        splicing_es_2 splicing_ad_2 splicing_aa_2 splicing_ada_2;
  by cell;
run;

proc print data=all_counts;
run;

%mend;

%counts(flag_immuno_gene);
%counts(flag_immunobase_diabetes_gene);

/* AI COUNTS:

             total_   ANNOT_   UNANNOT_                                                 ADA_
  cell        count    count     count    IR_count   ES_count   AD_count   AA_count    count

  cd4           490      176       136        178        156         41         51        71
  cd8           750      321       138        291        168         54         61        63
  cd19         1254      591       272        391        314        126        132       166
  cd4_8        1961      909       459        593        542        200        172       302
  cd4_19        224       98        55         71         52         18         24        24
  cd8_19        346      144        64        138         75         28         26        46
  cd4_8_19    20467    13824      3490       3153       3824       2146       2279      3730
  off        643481    10751    616121      16609     587021     130414     132902    114762


T1D COUNTS:

           total_   ANNOT_   UNANNOT_                                                ADA_
cell        count    count     count    IR_count   ES_count   AD_count   AA_count   count

cd4            95      39         22        34          21         12          5       16
cd8           177      66         37        74          38         11         19       17
cd19          306     134         82        90          95         37         39       51
cd4_8         491     264        111       116         127         62         32       79
cd4_19         41      18          3        20           6          2          6        1
cd8_19         68      25          7        36           9          3          5       14
cd4_8_19     4686    3087        829       770         803        493        510     1023
off        143122    2449     137171      3502      130195      28274      30301    30748

*/

/* Get DE and DD counts */

data splicing_clean_w_de_flags;
  set splicing_results_clean_w_flags;

  if flag_CD19_on=1 or flag_CD4_on=1 or flag_CD8_on=1 then flag_any_on=1;
  else flag_any_on=0;

  if flag_CD19_on=1 and flag_CD4_on=1 and flag_CD8_on=1 then flag_all_on=1;
  else flag_all_on=0;

  /* Tissue specificity */

  if flag_CD19_on=1 and flag_CD4_on=0 and flag_CD8_on=0 then flag_tissue_specific=1;
  else if flag_CD19_on=0 and flag_CD4_on=1 and flag_CD8_on=0 then flag_tissue_specific=1;
  else if flag_CD19_on=0 and flag_CD4_on=0 and flag_CD8_on=1 then flag_tissue_specificn=1;
  else flag_tissue_specific=0;

  /* Scenario 1 - only DE events expressed in all 3 tissues */
  
  if flag_all_on=1 and flag_anova_fdr_05=1 then flag_event_de=1;
  else if flag_all_on=1 and flag_anova_fdr_05=0 then flag_event_de=0;
  else if flag_all_on=0 then flag_event_de=.;

  /* Scenario 2 - only DE events expressed in all 3 tissues PLUS tissue-specific events */

  if flag_all_on=1 then do;
     if flag_anova_fdr_05=1 then flag_event_dd_1=1;
     else flag_event_dd_1=0;
     end;
  else if flag_all_on=0 then do;
      if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_event_dd_1=1;
      else if flag_cd4_on=0 and flag_cd8_on=1 and flag_cd19_on=0 then flag_event_dd_1=1;
      else if flag_cd4_on=0 and flag_cd8_on=0 and flag_cd19_on=1 then flag_event_dd_1=1;
      else flag_event_dd_1=.;
      end;

  /* Scenario 3 - only DE events expressed in all 3 tissues PLUS two-tissue events */


  if flag_all_on=1 then do;
     if flag_anova_fdr_05=1 then flag_event_dd_2=1;
     else flag_event_dd_2=0;
     end;
  else if flag_all_on=0 then do;
      if flag_cd4_on=1 or flag_cd8_on=1 or flag_cd19_on=1 then flag_event_dd_2=1;
      else flag_event_dd_2=.;
      end;

run;



/* Get DE and DD counts - diabetes genes */

%macro de_count(flag);

* total splicing events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* exon skipping events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_exonskip=1 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_alt_donor=1 and flag_alt_acceptor=0 and flag_intron_retention=0 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=0 and flag_alt_acceptor=1 and flag_intron_retention=0 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* alt donor and alt acceptor events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1  and flag_alt_donor=1 and flag_alt_acceptor=1 and flag_intron_retention=0 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* intron retention events;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=1 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

* annotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=1 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;


* unannotated junctions;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_junction_annotated=0 and flag_intron_retention=0 and &flag.=1;
   tables flag_event_de flag_event_dd_1 flag_event_dd_2;
run;

%mend;

%de_count(flag_immuno_gene);

/* AI DE COUNTS:
	DE	DD1	DD2
TOTAL	11775	14269	16800
ANNOT	8570	9658	10809
UNANNOT	2011	2557	3135
IR	1194	2054	2856
ES	2185	2823	3492
AD	1296	1517	1763
AA	1415	1659	1881
ADA	2289	2589	2961
*/

%de_count(flag_immunobase_diabetes_gene);


/* T1D DE COUNTS:
	DE	DD1	DD2
TOTAL	2466	3044	3644
ANNOT	1753	1992	2299
UNANNOT	438	579	700
IR	275	473	645
ES	424	578	720
AD	275	335	402
AA	297	360	403
ADA	584	668	762
*/


/* Redo splicing enrichments:
   autoimmune vs all genes
   T1D vs all autoimmune

   for:
   all events, IR, exonskip, unannotated */
   
   
/* Autoimmune among all detected */

proc freq data=splicing_clean_w_de_flags;
   tables flag_any_on*flag_immuno_gene / chisq ;
run;

/*
flag_any_on     flag_immuno_gene

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |5003158 | 643481 |5646639
         |  85.43 |  10.99 |  96.42
         |  88.60 |  11.40 |
         |  96.45 |  96.19 |
---------+--------+--------+
       1 | 184395 |  25492 | 209887
         |   3.15 |   0.44 |   3.58
         |  87.85 |  12.15 |
         |   3.55 |   3.81 |
---------+--------+--------+
Total     5187553   668973  5856526
            88.58    11.42   100.00

Statistics for Table of flag_any_on by flag_immuno_gene

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1    112.4317    <.0001
 Likelihood Ratio Chi-Square    1    110.5045    <.0001
 Continuity Adj. Chi-Square     1    112.3576    <.0001
 Mantel-Haenszel Chi-Square     1    112.4317    <.0001
 Phi Coefficient                       0.0044
 Contingency Coefficient               0.0044
 Cramer's V                            0.0044


                  Fisher's Exact Test
          -----------------------------------
          Cell (1,1) Frequency (F)    5003158
          Left-sided Pr <= F           1.0000
          Right-sided Pr >= F          <.0001

          Table Probability (P)        <.0001
          Two-sided Pr <= P            <.0001

                 Sample Size = 5856526
*/

/* T1D among autoimmune detected */

proc freq data=splicing_clean_w_de_flags;
   where flag_immuno_gene=1;
   tables flag_any_on*flag_immunobase_diabetes_gene / chisq ;
run;

/*
 flag_any_on
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 500359 | 143122 | 643481
          |  74.80 |  21.39 |  96.19
          |  77.76 |  22.24 |
          |  96.23 |  96.06 |
 ---------+--------+--------+
        1 |  19628 |   5864 |  25492
          |   2.93 |   0.88 |   3.81
          |  77.00 |  23.00 |
          |   3.77 |   3.94 |
 ---------+--------+--------+
 Total      519987   148986   668973
             77.73    22.27   100.00

Statistics for Table of flag_any_on by flag_immunobase_diabetes_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      8.2130    0.0042
       Likelihood Ratio Chi-Square    1      8.1524    0.0043
       Continuity Adj. Chi-Square     1      8.1691    0.0043
       Mantel-Haenszel Chi-Square     1      8.2130    0.0042
       Phi Coefficient                       0.0035
       Contingency Coefficient               0.0035
       Cramer's V                            0.0035


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)    500359
                 Left-sided Pr <= F          0.9979
                 Right-sided Pr >= F         0.0022

                 Table Probability (P)       0.0001
                 Two-sided Pr <= P           0.0043

                        Sample Size = 668973
*/


/* Tissue specificity */

data splicing_clean_w_de_flags2;
   set splicing_clean_w_de_flags;
   if flag_cd4_on=1 and flag_cd8_on=0 and flag_cd19_on=0 then flag_tissue_specific=1;
   else if flag_cd4_on=0 and flag_cd8_on=1 and flag_cd19_on=0 then flag_tissue_specific=1;
   else if flag_cd4_on=0 and flag_cd8_on=0 and flag_cd19_on=1 then flag_tissue_specific=1;
   else flag_tissue_specific=0;
run;

proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1;
   tables flag_tissue_specific*flag_immuno_gene / chisq;
run;

/*
 flag_tissue_specific
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 166054 |  22998 | 189052
          |  79.12 |  10.96 |  90.07
          |  87.84 |  12.16 |
          |  90.05 |  90.22 |
 ---------+--------+--------+
        1 |  18341 |   2494 |  20835
          |   8.74 |   1.19 |   9.93
          |  88.03 |  11.97 |
          |   9.95 |   9.78 |
 ---------+--------+--------+
 Total      184395    25492   209887
             87.85    12.15   100.00

Statistics for Table of flag_tissue_specific by flag_immuno_gene

     Statistic                     DF       Value      Prob
     ------------------------------------------------------
     Chi-Square                     1      0.6665    0.4143
     Likelihood Ratio Chi-Square    1      0.6689    0.4134
     Continuity Adj. Chi-Square     1      0.6484    0.4207
     Mantel-Haenszel Chi-Square     1      0.6665    0.4143
     Phi Coefficient                      -0.0018
     Contingency Coefficient               0.0018
     Cramer's V                           -0.0018


                      Fisher's Exact Test
               ----------------------------------
               Cell (1,1) Frequency (F)    166054
               Left-sided Pr <= F          0.2106
               Right-sided Pr >= F         0.7958

               Table Probability (P)       0.0064
               Two-sided Pr <= P           0.4211

                      Sample Size = 209887
*/


proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_tissue_specific*flag_immunobase_diabetes_gene / chisq;
run;

/*
 flag_tissue_specific
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  17712 |   5286 |  22998
          |  69.48 |  20.74 |  90.22
          |  77.02 |  22.98 |
          |  90.24 |  90.14 |
 ---------+--------+--------+
        1 |   1916 |    578 |   2494
          |   7.52 |   2.27 |   9.78
          |  76.82 |  23.18 |
          |   9.76 |   9.86 |
 ---------+--------+--------+
 Total       19628     5864    25492
             77.00    23.00   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      0.0464    0.8295
 Likelihood Ratio Chi-Square    1      0.0463    0.8297
 Continuity Adj. Chi-Square     1      0.0362    0.8491
 Mantel-Haenszel Chi-Square     1      0.0463    0.8295
 Phi Coefficient                       0.0013
 Contingency Coefficient               0.0013
 Cramer's V                            0.0013


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     17712
           Left-sided Pr <= F          0.5963
           Right-sided Pr >= F         0.4232

           Table Probability (P)       0.0195
           Two-sided Pr <= P           0.8217

                  Sample Size = 25492
*/


proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
run;

/*
Table of flag_event_dd_2 by flag_immuno_gene

    flag_event_dd_2
              flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  68143 |   8692 |  76835
             |  32.47 |   4.14 |  36.61
             |  88.69 |  11.31 |
             |  36.95 |  34.10 |
    ---------+--------+--------+
           1 | 116252 |  16800 | 133052
             |  55.39 |   8.00 |  63.39
             |  87.37 |  12.63 |
             |  63.05 |  65.90 |
    ---------+--------+--------+
    Total      184395    25492   209887
                87.85    12.15   100.00

tatistics for Table of flag_event_dd_2 by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     78.8248    <.0001
  Likelihood Ratio Chi-Square    1     79.5234    <.0001
  Continuity Adj. Chi-Square     1     78.7017    <.0001
  Mantel-Haenszel Chi-Square     1     78.8244    <.0001
  Phi Coefficient                       0.0194
  Contingency Coefficient               0.0194
  Cramer's V                            0.0194


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     68143
            Left-sided Pr <= F          1.0000
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 209887
*/


proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_immunobase_diabetes_gene / chisq;
run;

/*
  flag_event_dd_2
            flag_immunobase_diabetes_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   6472 |   2220 |   8692
           |  25.39 |   8.71 |  34.10
           |  74.46 |  25.54 |
           |  32.97 |  37.86 |
  ---------+--------+--------+
         1 |  13156 |   3644 |  16800
           |  51.61 |  14.29 |  65.90
           |  78.31 |  21.69 |
           |  67.03 |  62.14 |
  ---------+--------+--------+
  Total       19628     5864    25492
              77.00    23.00   100.00

             The SAS System       10:43 F

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     47.9447    <.0001
  Likelihood Ratio Chi-Square    1     47.4030    <.0001
  Continuity Adj. Chi-Square     1     47.7276    <.0001
  Mantel-Haenszel Chi-Square     1     47.9429    <.0001
  Phi Coefficient                      -0.0434
  Contingency Coefficient               0.0433
  Cramer's V                           -0.0434


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)      6472
            Left-sided Pr <= F          <.0001
            Right-sided Pr >= F         1.0000

            Table Probability (P)       <.0001
            Two-sided Pr <= P           <.0001

                   Sample Size = 25492

*/

/* Autoimmune among all DD2 */
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq ;
run;

/*
 Table of flag_event_dd_2 by flag_immuno_gene

     flag_event_dd_2
               flag_immuno_gene

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |  68143 |   8692 |  76835
              |  32.47 |   4.14 |  36.61
              |  88.69 |  11.31 |
              |  36.95 |  34.10 |
     ---------+--------+--------+
            1 | 116252 |  16800 | 133052
              |  55.39 |   8.00 |  63.39
              |  87.37 |  12.63 |
              |  63.05 |  65.90 |
     ---------+--------+--------+
     Total      184395    25492   209887
                 87.85    12.15   100.00

                The SAS System       10:43 Fri

              The FREQ Procedure
Statistics for Table of flag_event_dd_2 by flag_immuno_gen

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     78.8248    <.0001
   Likelihood Ratio Chi-Square    1     79.5234    <.0001
   Continuity Adj. Chi-Square     1     78.7017    <.0001
   Mantel-Haenszel Chi-Square     1     78.8244    <.0001
   Phi Coefficient                       0.0194
   Contingency Coefficient               0.0194
   Cramer's V                            0.0194


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)     68143
             Left-sided Pr <= F          1.0000
             Right-sided Pr >= F         <.0001

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 209887

*/

/* Diabetes among autoimmune DD2 */
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_immunobase_diabetes_gene / chisq ;
run;

/*
Table of flag_event_dd_2 by flag_immunobase_diabetes_gene

           flag_event_dd_2
                     flag_immunobase_diabetes_gene

           Frequency|
           Percent  |
           Row Pct  |
           Col Pct  |       0|       1|  Total
           ---------+--------+--------+
                  0 |   6472 |   2220 |   8692
                    |  25.39 |   8.71 |  34.10
                    |  74.46 |  25.54 |
                    |  32.97 |  37.86 |
           ---------+--------+--------+
                  1 |  13156 |   3644 |  16800
                    |  51.61 |  14.29 |  65.90
                    |  78.31 |  21.69 |
                    |  67.03 |  62.14 |
           ---------+--------+--------+
           Total       19628     5864    25492
                       77.00    23.00   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1     47.9447    <.0001
 Likelihood Ratio Chi-Square    1     47.4030    <.0001
 Continuity Adj. Chi-Square     1     47.7276    <.0001
 Mantel-Haenszel Chi-Square     1     47.9429    <.0001
 Phi Coefficient                      -0.0434
 Contingency Coefficient               0.0433
 Cramer's V                           -0.0434


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)      6472
           Left-sided Pr <= F          <.0001
           Right-sided Pr >= F         1.0000

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                  Sample Size = 25492
*/


/* Proportion of detected  ES events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1;
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
         0 |4494145 | 587021 |5081166
           |  87.65 |  11.45 |  99.10
           |  88.45 |  11.55 |
           |  99.10 |  99.13 |
  ---------+--------+--------+
         1 |  40950 |   5131 |  46081
           |   0.80 |   0.10 |   0.90
           |  88.87 |  11.13 |
           |   0.90 |   0.87 |
  ---------+--------+--------+
  Total     4535095   592152  5127247
              88.45    11.55   100.00

 Statistics for Table of flag_any_on by flag_immuno_gene

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1      7.8161    0.0052
  Likelihood Ratio Chi-Square    1      7.8982    0.0049
  Continuity Adj. Chi-Square     1      7.7752    0.0053
  Mantel-Haenszel Chi-Square     1      7.8161    0.0052
  Phi Coefficient                      -0.0012
  Contingency Coefficient               0.0012
  Cramer's V                           -0.0012


                   Fisher's Exact Test
           -----------------------------------
           Cell (1,1) Frequency (F)    4494145
           Left-sided Pr <= F           0.0025
           Right-sided Pr >= F          0.9976

           Table Probability (P)        0.0001
           Two-sided Pr <= P            0.0051

                  Sample Size = 5127247
*/


/* Proportion of detected autoimmune ES events that are diabetes events */


proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_immuno_gene=1;
   tables flag_any_on*flag_immunobase_diabetes_gene / chisq;
run;


/*
 Table of flag_any_on by flag_immunobase_diabetes_gene

          flag_any_on
                    flag_immunobase_diabetes_gene

          Frequency|
          Percent  |
          Row Pct  |
          Col Pct  |       0|       1|  Total
          ---------+--------+--------+
                 0 | 456826 | 130195 | 587021
                   |  77.15 |  21.99 |  99.13
                   |  77.82 |  22.18 |
                   |  99.13 |  99.16 |
          ---------+--------+--------+
                 1 |   4032 |   1099 |   5131
                   |   0.68 |   0.19 |   0.87
                   |  78.58 |  21.42 |
                   |   0.87 |   0.84 |
          ---------+--------+--------+
          Total      460858   131294   592152
                      77.83    22.17   100.00

                     The SAS System       10:43 Friday,

Statistics for Table of flag_any_on by flag_immunobase_diabetes_gene

       Statistic                     DF       Value      Prob
       ------------------------------------------------------
       Chi-Square                     1      1.7030    0.1919
       Likelihood Ratio Chi-Square    1      1.7170    0.1901
       Continuity Adj. Chi-Square     1      1.6593    0.1977
       Mantel-Haenszel Chi-Square     1      1.7030    0.1919
       Phi Coefficient                      -0.0017
       Contingency Coefficient               0.0017
       Cramer's V                           -0.0017


                        Fisher's Exact Test
                 ----------------------------------
                 Cell (1,1) Frequency (F)    456826
                 Left-sided Pr <= F          0.0985
                 Right-sided Pr >= F         0.9073

                 Table Probability (P)       0.0058
                 Two-sided Pr <= P           0.1938

                        Sample Size = 592152
*/


/* Proportion of  DE ES events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
run;


/*
Table of flag_event_dd_2 by flag_immuno_gene

    flag_event_dd_2
              flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |  14771 |   1639 |  16410
             |  32.05 |   3.56 |  35.61
             |  90.01 |   9.99 |
             |  36.07 |  31.94 |
    ---------+--------+--------+
           1 |  26179 |   3492 |  29671
             |  56.81 |   7.58 |  64.39
             |  88.23 |  11.77 |
             |  63.93 |  68.06 |
    ---------+--------+--------+
    Total       40950     5131    46081
                88.87    11.13   100.00

 Statistics for Table of flag_event_dd_2 by flag_immuno_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     33.8812    <.0001
    Likelihood Ratio Chi-Square    1     34.3782    <.0001
    Continuity Adj. Chi-Square     1     33.7014    <.0001
    Mantel-Haenszel Chi-Square     1     33.8805    <.0001
    Phi Coefficient                       0.0271
    Contingency Coefficient               0.0271
    Cramer's V                            0.0271


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     14771
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 46081
*/

/* Proportion of autoimmune DE ES events that are T1D events */

proc freq data=splicing_clean_w_de_flags;
   where flag_exonskip=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_immunobase_diabetes_gene / chisq;
run;


/*
 flag_event_dd_2
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

            The SAS System       10:43 Fri


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

                  Sample Size = 5131

*/


/* Proportion of detected  IR events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1;
   tables flag_any_on*flag_immuno_gene / chisq;
run;

/*
 flag_any_on     flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 167679 |  16609 | 184288
          |  74.46 |   7.38 |  81.84
          |  90.99 |   9.01 |
          |  82.29 |  77.53 |
 ---------+--------+--------+
        1 |  36077 |   4815 |  40892
          |  16.02 |   2.14 |  18.16
          |  88.23 |  11.77 |
          |  17.71 |  22.47 |
 ---------+--------+--------+
 Total      203756    21424   225180
             90.49     9.51   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1    296.6374    <.0001
 Likelihood Ratio Chi-Square    1    282.1158    <.0001
 Continuity Adj. Chi-Square     1    296.3166    <.0001
 Mantel-Haenszel Chi-Square     1    296.6361    <.0001
 Phi Coefficient                       0.0363
 Contingency Coefficient               0.0363
 Cramer's V                            0.0363


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)    167679
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                  Sample Size = 225180

*/

/* Proportion of detected autoimmune IR events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_immuno_gene=1;
   tables flag_any_on*flag_immunobase_diabetes_gene / chisq;
run;

/*
 flag_any_on
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  13107 |   3502 |  16609
          |  61.18 |  16.35 |  77.53
          |  78.92 |  21.08 |
          |  78.10 |  75.44 |
 ---------+--------+--------+
        1 |   3675 |   1140 |   4815
          |  17.15 |   5.32 |  22.47
          |  76.32 |  23.68 |
          |  21.90 |  24.56 |
 ---------+--------+--------+
 Total       16782     4642    21424
             78.33    21.67   100.00

  Statistic                     DF       Value      Prob
  ------------------------------------------------------
  Chi-Square                     1     14.7654    0.0001
  Likelihood Ratio Chi-Square    1     14.5438    0.0001
  Continuity Adj. Chi-Square     1     14.6132    0.0001
  Mantel-Haenszel Chi-Square     1     14.7647    0.0001
  Phi Coefficient                       0.0263
  Contingency Coefficient               0.0262
  Cramer's V                            0.0263


                   Fisher's Exact Test
            ----------------------------------
            Cell (1,1) Frequency (F)     13107
            Left-sided Pr <= F          0.9999
            Right-sided Pr >= F         <.0001

            Table Probability (P)       <.0001
            Two-sided Pr <= P           0.0001

                   Sample Size = 21424

*/

/* Proportion of  DE IR events that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
run;

/*
 flag_event_dd_2
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  15297 |   1959 |  17256
          |  37.41 |   4.79 |  42.20
          |  88.65 |  11.35 |
          |  42.40 |  40.69 |
 ---------+--------+--------+
        1 |  20780 |   2856 |  23636
          |  50.82 |   6.98 |  57.80
          |  87.92 |  12.08 |
          |  57.60 |  59.31 |
 ---------+--------+--------+
 Total       36077     4815    40892
             88.23    11.77   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      5.1262    0.0236
 Likelihood Ratio Chi-Square    1      5.1413    0.0234
 Continuity Adj. Chi-Square     1      5.0561    0.0245
 Mantel-Haenszel Chi-Square     1      5.1260    0.0236
 Phi Coefficient                       0.0112
 Contingency Coefficient               0.0112
 Cramer's V                            0.0112


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     15297
           Left-sided Pr <= F          0.9888
           Right-sided Pr >= F         0.0122

           Table Probability (P)       0.0010
           Two-sided Pr <= P           0.0243

                  Sample Size = 40892
*/

/* Proportion of autoimmune DE IR events that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=1 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_immunobase_diabetes_gene / chisq;
run;


/*
 flag_event_dd_2
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   1464 |    495 |   1959
          |  30.40 |  10.28 |  40.69
          |  74.73 |  25.27 |
          |  39.84 |  43.42 |
 ---------+--------+--------+
        1 |   2211 |    645 |   2856
          |  45.92 |  13.40 |  59.31
          |  77.42 |  22.58 |
          |  60.16 |  56.58 |
 ---------+--------+--------+
 Total        3675     1140     4815
             76.32    23.68   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      4.6321    0.0314
 Likelihood Ratio Chi-Square    1      4.6117    0.0318
 Continuity Adj. Chi-Square     1      4.4848    0.0342
 Mantel-Haenszel Chi-Square     1      4.6311    0.0314
 Phi Coefficient                      -0.0310
 Contingency Coefficient               0.0310
 Cramer's V                           -0.0310


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)      1464
           Left-sided Pr <= F          0.0173
           Right-sided Pr >= F         0.9855

           Table Probability (P)       0.0027
           Two-sided Pr <= P           0.0324

                   Sample Size = 4815

*/

/* Among detected junctions in all genes, what is the proportion of unannotated from autoimmune genes ?*/

proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=0;
   tables flag_junction_annotated*flag_immuno_gene / chisq ;
run;

/*

Table of flag_junction_annotated by flag_immuno_gene

        flag_junction_annotated
                  flag_immuno_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |  34103 |   4614 |  38717
                 |  20.18 |   2.73 |  22.91
                 |  88.08 |  11.92 |
                 |  22.99 |  22.31 |
        ---------+--------+--------+
               1 | 114215 |  16063 | 130278
                 |  67.58 |   9.51 |  77.09
                 |  87.67 |  12.33 |
                 |  77.01 |  77.69 |
        ---------+--------+--------+
        Total      148318    20677   168995
                    87.76    12.24   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      4.7304    0.0296
 Likelihood Ratio Chi-Square    1      4.7556    0.0292
 Continuity Adj. Chi-Square     1      4.6921    0.0303
 Mantel-Haenszel Chi-Square     1      4.7304    0.0296
 Phi Coefficient                       0.0053
 Contingency Coefficient               0.0053
 Cramer's V                            0.0053


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     34103
           Left-sided Pr <= F          0.9857
           Right-sided Pr >= F         0.0150

           Table Probability (P)       0.0007
           Two-sided Pr <= P           0.0298

                  Sample Size = 168995
*/

/* Among detected junctions in autoimmune genes, what is the proportion of unannotated from diabetes genes ?*/
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_junction_annotated*flag_immunobase_diabetes_gene / chisq ;
run;


/*
 flag_junction_annotated
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   3523 |   1091 |   4614
          |  17.04 |   5.28 |  22.31
          |  76.35 |  23.65 |
          |  22.08 |  23.09 |
 ---------+--------+--------+
        1 |  12430 |   3633 |  16063
          |  60.12 |  17.57 |  77.69
          |  77.38 |  22.62 |
          |  77.92 |  76.91 |
 ---------+--------+--------+
 Total       15953     4724    20677
             77.15    22.85   100.00

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      2.1499    0.1426
Likelihood Ratio Chi-Square    1      2.1376    0.1437
Continuity Adj. Chi-Square     1      2.0920    0.1481
Mantel-Haenszel Chi-Square     1      2.1498    0.1426
Phi Coefficient                      -0.0102
Contingency Coefficient               0.0102
Cramer's V                           -0.0102


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)      3523
          Left-sided Pr <= F          0.0743
          Right-sided Pr >= F         0.9311

          Table Probability (P)       0.0054
          Two-sided Pr <= P           0.1465

                 Sample Size = 20677
*/

/* Among detected events in  genes, what is the proportion of IR events from autoimmune genes ?*/
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_intron_retention*flag_immuno_gene / chisq ;
run;

/*
 Table of flag_intron_retention by flag_immuno_gene

        flag_intron_retention
                  flag_immuno_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 | 148318 |  20677 | 168995
                 |  70.67 |   9.85 |  80.52
                 |  87.76 |  12.24 |
                 |  80.43 |  81.11 |
        ---------+--------+--------+
               1 |  36077 |   4815 |  40892
                 |  17.19 |   2.29 |  19.48
                 |  88.23 |  11.77 |
                 |  19.57 |  18.89 |
        ---------+--------+--------+
        Total      184395    25492   209887
                    87.85    12.15   100.00

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      6.5393    0.0106
Likelihood Ratio Chi-Square    1      6.5835    0.0103
Continuity Adj. Chi-Square     1      6.4962    0.0108
Mantel-Haenszel Chi-Square     1      6.5392    0.0106
Phi Coefficient                      -0.0056
Contingency Coefficient               0.0056
Cramer's V                           -0.0056


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)    148318
          Left-sided Pr <= F          0.0053
          Right-sided Pr >= F         0.9950

          Table Probability (P)       0.0003
          Two-sided Pr <= P           0.0106

                 Sample Size = 209887
*/

/* Among detected events in autoimmune genes, what is the proportion of IR events from diabetes genes ?*/
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_intron_retention*flag_immunobase_diabetes_gene / chisq ;
run;


/*
 flag_intron_retention
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  15953 |   4724 |  20677
          |  62.58 |  18.53 |  81.11
          |  77.15 |  22.85 |
          |  81.28 |  80.56 |
 ---------+--------+--------+
        1 |   3675 |   1140 |   4815
          |  14.42 |   4.47 |  18.89
          |  76.32 |  23.68 |
          |  18.72 |  19.44 |
 ---------+--------+--------+
 Total       19628     5864    25492
             77.00    23.00   100.00

            The SAS System       10:43 Frid

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      1.5168    0.2181
Likelihood Ratio Chi-Square    1      1.5089    0.2193
Continuity Adj. Chi-Square     1      1.4703    0.2253
Mantel-Haenszel Chi-Square     1      1.5167    0.2181
Phi Coefficient                       0.0077
Contingency Coefficient               0.0077
Cramer's V                            0.0077


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)     15953
          Left-sided Pr <= F          0.8942
          Right-sided Pr >= F         0.1128

          Table Probability (P)       0.0071
          Two-sided Pr <= P           0.2237

                 Sample Size = 25492

*/

/* Proportion of DE unannotated junctions among all genes that are autoimmune events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_any_on=1;
   tables flag_event_dd_2*flag_immuno_gene / chisq;
run;

/*
  flag_event_dd_2
            flag_immuno_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  12280 |   1479 |  13759
           |  31.72 |   3.82 |  35.54
           |  89.25 |  10.75 |
           |  36.01 |  32.05 |
  ---------+--------+--------+
         1 |  21823 |   3135 |  24958
           |  56.37 |   8.10 |  64.46
           |  87.44 |  12.56 |
           |  63.99 |  67.95 |
  ---------+--------+--------+
  Total       34103     4614    38717
              88.08    11.92   100.00

 Statistics for Table of flag_event_dd_2 by flag_immuno_gene

    Statistic                     DF       Value      Prob
    ------------------------------------------------------
    Chi-Square                     1     27.7355    <.0001
    Likelihood Ratio Chi-Square    1     28.1181    <.0001
    Continuity Adj. Chi-Square     1     27.5632    <.0001
    Mantel-Haenszel Chi-Square     1     27.7348    <.0001
    Phi Coefficient                       0.0268
    Contingency Coefficient               0.0268
    Cramer's V                            0.0268


                     Fisher's Exact Test
              ----------------------------------
              Cell (1,1) Frequency (F)     12280
              Left-sided Pr <= F          1.0000
              Right-sided Pr >= F         <.0001

              Table Probability (P)       <.0001
              Two-sided Pr <= P           <.0001

                     Sample Size = 38717
*/


/* Proportion of DE unannotated junctions among autoimmune genes that are diabetes events */

proc freq data=splicing_clean_w_de_flags;
   where flag_intron_retention=0 and flag_junction_annotated=0 and flag_any_on=1 and flag_immuno_gene=1;
   tables flag_event_dd_2*flag_immunobase_diabetes_gene / chisq;
run;


/*

flag_event_dd_2
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

                   Sample Size = 4614


*/


/* Proportion of exon skipping of all detected events */


* AUTOIMMUNE vs ALL;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_exonskip*flag_immuno_gene / chisq;
run;

/*
 Table of flag_exonskip by flag_immuno_gene

    flag_exonskip     flag_immuno_gene

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 | 143445 |  20361 | 163806
             |  68.34 |   9.70 |  78.04
             |  87.57 |  12.43 |
             |  77.79 |  79.87 |
    ---------+--------+--------+
           1 |  40950 |   5131 |  46081
             |  19.51 |   2.44 |  21.96
             |  88.87 |  11.13 |
             |  22.21 |  20.13 |
    ---------+--------+--------+
    Total      184395    25492   209887
                87.85    12.15   100.00

               The SAS System       10:43 Friday,


 Statistics for Table of flag_exonskip by flag_immuno_gene

   Statistic                     DF       Value      Prob
   ------------------------------------------------------
   Chi-Square                     1     56.5408    <.0001
   Likelihood Ratio Chi-Square    1     57.5605    <.0001
   Continuity Adj. Chi-Square     1     56.4195    <.0001
   Mantel-Haenszel Chi-Square     1     56.5406    <.0001
   Phi Coefficient                      -0.0164
   Contingency Coefficient               0.0164
   Cramer's V                           -0.0164


                    Fisher's Exact Test
             ----------------------------------
             Cell (1,1) Frequency (F)    143445
             Left-sided Pr <= F          <.0001
             Right-sided Pr >= F         1.0000

             Table Probability (P)       <.0001
             Two-sided Pr <= P           <.0001

                    Sample Size = 209887
*/


* T1D vs AUTOIMMUNE ;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_exonskip*flag_immunobase_diabetes_gene / chisq;
run;

/*

  flag_exonskip
            flag_immunobase_diabetes_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  15596 |   4765 |  20361
           |  61.18 |  18.69 |  79.87
           |  76.60 |  23.40 |
           |  79.46 |  81.26 |
  ---------+--------+--------+
         1 |   4032 |   1099 |   5131
           |  15.82 |   4.31 |  20.13
           |  78.58 |  21.42 |
           |  20.54 |  18.74 |
  ---------+--------+--------+
  Total       19628     5864    25492
              77.00    23.00   100.00

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      9.1057    0.0025
Likelihood Ratio Chi-Square    1      9.2201    0.0024
Continuity Adj. Chi-Square     1      8.9940    0.0027
Mantel-Haenszel Chi-Square     1      9.1053    0.0025
Phi Coefficient                      -0.0189
Contingency Coefficient               0.0189
Cramer's V                           -0.0189


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)     15596
          Left-sided Pr <= F          0.0013
          Right-sided Pr >= F         0.9989

          Table Probability (P)       0.0002
          Two-sided Pr <= P           0.0025

                 Sample Size = 25492

*/

/* Proportion of IR events of all detected events */

* AUTOIMMUNE vs ALL;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_intron_retention*flag_immuno_gene / chisq;
run;

/*
 Table of flag_intron_retention by flag_immuno_gene

        flag_intron_retention
                  flag_immuno_gene

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 | 148318 |  20677 | 168995
                 |  70.67 |   9.85 |  80.52
                 |  87.76 |  12.24 |
                 |  80.43 |  81.11 |
        ---------+--------+--------+
               1 |  36077 |   4815 |  40892
                 |  17.19 |   2.29 |  19.48
                 |  88.23 |  11.77 |
                 |  19.57 |  18.89 |
        ---------+--------+--------+
        Total      184395    25492   209887
                    87.85    12.15   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      6.5393    0.0106
 Likelihood Ratio Chi-Square    1      6.5835    0.0103
 Continuity Adj. Chi-Square     1      6.4962    0.0108
 Mantel-Haenszel Chi-Square     1      6.5392    0.0106
 Phi Coefficient                      -0.0056
 Contingency Coefficient               0.0056
 Cramer's V                           -0.0056


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)    148318
           Left-sided Pr <= F          0.0053
           Right-sided Pr >= F         0.9950

           Table Probability (P)       0.0003
           Two-sided Pr <= P           0.0106

                  Sample Size = 209887
*/

* T1D vs AUTOIMMUNE ;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_intron_retention*flag_immunobase_diabetes_gene / chisq;
run;

/*
flag_intron_retention
          flag_immunobase_diabetes_gene

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  15953 |   4724 |  20677
         |  62.58 |  18.53 |  81.11
         |  77.15 |  22.85 |
         |  81.28 |  80.56 |
---------+--------+--------+
       1 |   3675 |   1140 |   4815
         |  14.42 |   4.47 |  18.89
         |  76.32 |  23.68 |
         |  18.72 |  19.44 |
---------+--------+--------+
Total       19628     5864    25492
            77.00    23.00   100.00

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      1.5168    0.2181
Likelihood Ratio Chi-Square    1      1.5089    0.2193
Continuity Adj. Chi-Square     1      1.4703    0.2253
Mantel-Haenszel Chi-Square     1      1.5167    0.2181
Phi Coefficient                       0.0077
Contingency Coefficient               0.0077
Cramer's V                            0.0077


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)     15953
          Left-sided Pr <= F          0.8942
          Right-sided Pr >= F         0.1128

          Table Probability (P)       0.0071
          Two-sided Pr <= P           0.2237

                 Sample Size = 25492
*/


/* Proportion of Unannotated events of all detected events */


data splicing_clean_w_de_flags3;
   set splicing_clean_w_de_flags;
   if flag_junction_annotated=0 and flag_intron_retention=0 then flaG_junction_unannotated=1;
   else flag_junction_unannotated=0;
run;

* AUTOIMMUNE vs ALL;
proc freq data=splicing_clean_w_de_flags3;
   where flag_any_on=1;
   tables flag_junction_unannotated*flag_immuno_gene / chisq ;
run;

/*
   flaG_junction_unannotated
             flag_immuno_gene

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 | 150292 |  20878 | 171170
            |  71.61 |   9.95 |  81.55
            |  87.80 |  12.20 |
            |  81.51 |  81.90 |
   ---------+--------+--------+
          1 |  34103 |   4614 |  38717
            |  16.25 |   2.20 |  18.45
            |  88.08 |  11.92 |
            |  18.49 |  18.10 |
   ---------+--------+--------+
   Total      184395    25492   209887
               87.85    12.15   100.00

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      2.3197    0.1277
Likelihood Ratio Chi-Square    1      2.3295    0.1269
Continuity Adj. Chi-Square     1      2.2935    0.1299
Mantel-Haenszel Chi-Square     1      2.3197    0.1277
Phi Coefficient                      -0.0033
Contingency Coefficient               0.0033
Cramer's V                           -0.0033


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)    150292
          Left-sided Pr <= F          0.0647
          Right-sided Pr >= F         0.9374

          Table Probability (P)       0.0022
          Two-sided Pr <= P           0.1295

                 Sample Size = 209887


*/

* T1D vs AUTOIMMUNE ;
proc freq data=splicing_clean_w_de_flags3;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_junction_unannotated*flag_immunobase_diabetes_gene / chisq ;
run;

/*
 flaG_junction_unannotated
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  16105 |   4773 |  20878
          |  63.18 |  18.72 |  81.90
          |  77.14 |  22.86 |
          |  82.05 |  81.39 |
 ---------+--------+--------+
        1 |   3523 |   1091 |   4614
          |  13.82 |   4.28 |  18.10
          |  76.35 |  23.65 |
          |  17.95 |  18.61 |
 ---------+--------+--------+
 Total       19628     5864    25492
             77.00    23.00   100.00


 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      1.3115    0.2521
 Likelihood Ratio Chi-Square    1      1.3050    0.2533
 Continuity Adj. Chi-Square     1      1.2676    0.2602
 Mantel-Haenszel Chi-Square     1      1.3115    0.2521
 Phi Coefficient                       0.0072
 Contingency Coefficient               0.0072
 Cramer's V                            0.0072


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     16105
           Left-sided Pr <= F          0.8777
           Right-sided Pr >= F         0.1302

           Table Probability (P)       0.0080
           Two-sided Pr <= P           0.2542

                  Sample Size = 25492

*/

/* IR+Unannot vs Annotated junctions */

* AUTOIMMUNE vs ALL;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1;
   tables flag_junction_annotated*flag_immuno_gene / chisq ;
run;

/*
  flag_junction_annotated
            flag_immuno_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  70180 |   9429 |  79609
           |  33.44 |   4.49 |  37.93
           |  88.16 |  11.84 |
           |  38.06 |  36.99 |
  ---------+--------+--------+
         1 | 114215 |  16063 | 130278
           |  54.42 |   7.65 |  62.07
           |  87.67 |  12.33 |
           |  61.94 |  63.01 |
  ---------+--------+--------+
  Total      184395    25492   209887
              87.85    12.15   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1     10.9222    0.0010
 Likelihood Ratio Chi-Square    1     10.9533    0.0009
 Continuity Adj. Chi-Square     1     10.8767    0.0010
 Mantel-Haenszel Chi-Square     1     10.9222    0.0010
 Phi Coefficient                       0.0072
 Contingency Coefficient               0.0072
 Cramer's V                            0.0072


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     70180
           Left-sided Pr <= F          0.9995
           Right-sided Pr >= F         0.0005

           Table Probability (P)       <.0001
           Two-sided Pr <= P           0.0009

                  Sample Size = 209887



*/


* T1D vs AUTOIMMUNE;
proc freq data=splicing_clean_w_de_flags;
   where flag_any_on=1 and flag_immuno_gene=1;
   tables flag_junction_annotated*flag_immunobase_diabetes_gene / chisq ;
run;

/*
 flag_junction_annotated
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   7198 |   2231 |   9429
          |  28.24 |   8.75 |  36.99
          |  76.34 |  23.66 |
          |  36.67 |  38.05 |
 ---------+--------+--------+
        1 |  12430 |   3633 |  16063
          |  48.76 |  14.25 |  63.01
          |  77.38 |  22.62 |
          |  63.33 |  61.95 |
 ---------+--------+--------+
 Total       19628     5864    25492
             77.00    23.00   100.00

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      3.6551    0.0559
Likelihood Ratio Chi-Square    1      3.6453    0.0562
Continuity Adj. Chi-Square     1      3.5964    0.0579
Mantel-Haenszel Chi-Square     1      3.6550    0.0559
Phi Coefficient                      -0.0120
Contingency Coefficient               0.0120
Cramer's V                           -0.0120


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)      7198
          Left-sided Pr <= F          0.0291
          Right-sided Pr >= F         0.9729

          Table Probability (P)       0.0020
          Two-sided Pr <= P           0.0560

                 Sample Size = 25492

*/


/* Tissue spec - ES, IR, Unannot */


*ES - AI vs other;
proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_exonskip=1;
   tables flag_tissue_specific*flag_immuno_gene / chisq ;
run;

/*
flag_tissue_specific
          flag_immuno_gene

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  36155 |   4493 |  40648
         |  78.46 |   9.75 |  88.21
         |  88.95 |  11.05 |
         |  88.29 |  87.57 |
---------+--------+--------+
       1 |   4795 |    638 |   5433
         |  10.41 |   1.38 |  11.79
         |  88.26 |  11.74 |
         |  11.71 |  12.43 |
---------+--------+--------+
Total       40950     5131    46081
            88.87    11.13   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      2.3034    0.1291
 Likelihood Ratio Chi-Square    1      2.2724    0.1317
 Continuity Adj. Chi-Square     1      2.2342    0.1350
 Mantel-Haenszel Chi-Square     1      2.3033    0.1291
 Phi Coefficient                       0.0071
 Contingency Coefficient               0.0071
 Cramer's V                            0.0071


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     36155
           Left-sided Pr <= F          0.9376
           Right-sided Pr >= F         0.0682

           Table Probability (P)       0.0058
           Two-sided Pr <= P           0.1297

                  Sample Size = 46081
*/

*ES - T1D vs Other AI;
proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_exonskip=1 and flag_immuno_gene=1;
   tables flag_tissue_specific*flag_immunobase_diabetes_gene / chisq ;
run;

/*

 flag_tissue_specific
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

                  Sample Size = 5131

*/

*IR - AI vs other;
proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_intron_retention=1;
   tables flag_tissue_specific*flag_immuno_gene / chisq ;
run;

/*
flag_tissue_specific
          flag_immuno_gene

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |  30310 |   3955 |  34265
         |  74.12 |   9.67 |  83.79
         |  88.46 |  11.54 |
         |  84.01 |  82.14 |
---------+--------+--------+
       1 |   5767 |    860 |   6627
         |  14.10 |   2.10 |  16.21
         |  87.02 |  12.98 |
         |  15.99 |  17.86 |
---------+--------+--------+
Total       36077     4815    40892
            88.23    11.77   100.00

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1     11.0047    0.0009
 Likelihood Ratio Chi-Square    1     10.7563    0.0010
 Continuity Adj. Chi-Square     1     10.8670    0.0010
 Mantel-Haenszel Chi-Square     1     11.0044    0.0009
 Phi Coefficient                       0.0164
 Contingency Coefficient               0.0164
 Cramer's V                            0.0164


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)     30310
           Left-sided Pr <= F          0.9995
           Right-sided Pr >= F         0.0006

           Table Probability (P)       <.0001
           Two-sided Pr <= P           0.0010

                  Sample Size = 40892

*/

*IR - T1D vs Other AI;
proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_intron_retention=1 and flag_immuno_gene=1;
   tables flag_tissue_specific*flag_immunobase_diabetes_gene / chisq ;
run;

/*
           flag_immunobase_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   3013 |    942 |   3955
          |  62.58 |  19.56 |  82.14
          |  76.18 |  23.82 |
          |  81.99 |  82.63 |
 ---------+--------+--------+
        1 |    662 |    198 |    860
          |  13.75 |   4.11 |  17.86
          |  76.98 |  23.02 |
          |  18.01 |  17.37 |
 ---------+--------+--------+
 Total        3675     1140     4815
             76.32    23.68   100.00


 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      0.2469    0.6193
 Likelihood Ratio Chi-Square    1      0.2481    0.6184
 Continuity Adj. Chi-Square     1      0.2049    0.6508
 Mantel-Haenszel Chi-Square     1      0.2468    0.6193
 Phi Coefficient                      -0.0072
 Contingency Coefficient               0.0072
 Cramer's V                           -0.0072


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)      3013
           Left-sided Pr <= F          0.3269
           Right-sided Pr >= F         0.7046

           Table Probability (P)       0.0314
           Two-sided Pr <= P           0.6581

                   Sample Size = 4815
*/

*Unannot - AI vs other;
proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_junction_annotated=0 and flag_intron_retention=0;
   tables flag_tissue_specific*flag_immuno_gene / chisq ;
run;

/*
 flag_tissue_specific
           flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  30015 |   4068 |  34083
          |  77.52 |  10.51 |  88.03
          |  88.06 |  11.94 |
          |  88.01 |  88.17 |
 ---------+--------+--------+
        1 |   4088 |    546 |   4634
          |  10.56 |   1.41 |  11.97
          |  88.22 |  11.78 |
          |  11.99 |  11.83 |
 ---------+--------+--------+
 Total       34103     4614    38717
             88.08    11.92   100.00

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1      0.0911    0.7628
Likelihood Ratio Chi-Square    1      0.0913    0.7625
Continuity Adj. Chi-Square     1      0.0771    0.7813
Mantel-Haenszel Chi-Square     1      0.0911    0.7628
Phi Coefficient                      -0.0015
Contingency Coefficient               0.0015
Cramer's V                           -0.0015


                 Fisher's Exact Test
          ----------------------------------
          Cell (1,1) Frequency (F)     30015
          Left-sided Pr <= F          0.3923
          Right-sided Pr >= F         0.6262

          Table Probability (P)       0.0185
          Two-sided Pr <= P           0.7904

                 Sample Size = 38717





*/

*Unannot - T1D vs Other AI;
proc freq data=splicing_clean_w_de_flags2;
   where flag_any_on=1 and flag_junction_annotated=0 and flag_intron_retention=0 and flag_immuno_gene=1;
   tables flag_tissue_specific*flag_immunobase_diabetes_gene / chisq ;
run;


/*
 flag_tissue_specific
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

            The SAS System       10:43 F

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

                  Sample Size = 4614
*/

quit;


