/* Get counts for DE exons */


data results_by_fusion;
   set con.results_by_fusion_final2;
run;


/* Need to count:
#fusions expressed per Venn group
#genes expressed per Venn group (same as eQTL counts) --> only single-gene fusions
#fusions DE (each comparison, then general ANOVA
#genes DE (each comparison, then general ANOVA --> only single-gene fusions

*/


/**** FUSIONS EXPRESSED PER CELL TYPE GROUP *****/

data fusions_expressed;
   set results_by_fusion;
   length cell_group $12.;
   if flag_CD4_on=1 and flag_CD8_on=0 and flag_CD19_on=0 then cell_group='CD4';
   else if flag_CD4_on=0 and flag_CD8_on=1 and flag_CD19_on=0 then cell_group='CD8';
   else if flag_CD4_on=0 and flag_CD8_on=0 and flag_CD19_on=1 then cell_group='CD19';
   else if flag_CD4_on=1 and flag_CD8_on=1 and flag_CD19_on=0 then cell_group='CD4_CD8';
   else if flag_CD4_on=1 and flag_CD8_on=0 and flag_CD19_on=1 then cell_group='CD4_CD19';
   else if flag_CD4_on=0 and flag_CD8_on=1 and flag_CD19_on=1 then cell_group='CD8_CD19';
   else if flag_CD4_on=1 and flag_CD8_on=1 and flag_CD19_on=1 then cell_group='CD4_CD8_CD19';
   else cell_group='None';
   keep fusion_id gene_cat flag_CD4_on flag_CD8_on flag_CD19_on cell_group;
run;

proc freq data=fusions_expressed;
   tables cell_group;
run;


/****** GENES EXPRESSED PER CELL TYPE GROUP *******/

/* ALL GENES */
/* Split genes on cell type */

data immunogene_flags;
   set con.immunogene_flags;
run;

proc sort data=immunogene_flags;
   by gene_id;
run;

data results_by_fusion2;
   set results_by_fusion;
   rename gene_cat=gene_id;
run;

proc sort data=results_by_fusion2;
   by gene_id;
run;

data results_by_fusion_w_flags;
   merge results_by_fusion2 (in=in1) immunogene_flags (in=in2);
   by gene_id;
   if in1;
run;

data CD4_genes;
   set results_by_fusion_w_flags;
   if flag_CD4_on=1 and flag_multigene=0 and flag_pseudogene=0;
   keep gene_id;
run;

data CD8_genes;
   set results_by_fusion_w_flags;
   if flag_CD8_on=1 and flag_multigene=0 and flag_pseudogene=0;
   keep gene_id;
run;

data CD19_genes;
   set results_by_fusion_w_flags;
   if flag_CD19_on=1 and flag_multigene=0 and flag_pseudogene=0;
   keep gene_id;
run;

data all_genes;
   set results_by_fusion_w_flags;
   if flag_multigene=0 and flag_pseudogene=0;
   keep gene_id;
run;

proc sort data=CD4_genes nodup;
  by gene_id;
run;

proc sort data=CD8_genes nodup;
  by gene_id;
run;

proc sort data=CD19_genes nodup;
  by gene_id;
run;

proc sort data=all_genes nodup;
  by gene_id;
run;

data expressed_genes_by_group;
   merge CD4_genes (in=in1) CD8_genes (in=in2) CD19_genes (in=in3) all_genes (in=in4);
   by gene_id;
   length cell_group $12.;
   if in1 and in2 and in3 then cell_group='CD4_CD8_CD19';
 
   else if in1 and not in2 and not in3 then cell_group='CD4';
   else if not in1 and in2 and not in3 then cell_group='CD8';
   else if not in1 and not in2 and in3 then cell_group='CD19';

   else if in1 and in2 and not in3 then cell_group='CD4_CD8';
   else if in1 and not in2 and in3 then cell_group='CD4_CD19';
   else if not in1 and in2 and in3 then cell_group='CD8_CD19';

   else cell_group='OFF';
run;

proc freq data=expressed_genes_by_group;
   tables cell_group;
run;

/***** NOW DO FOR T1D and IMMUNO GENES ****/

/* FUSIONS */
/* T1D genes */


data t1d_fusions_expressed;
   set results_by_fusion_w_flags;
   length cell_group $12.;
   if flag_CD4_on=1 and flag_CD8_on=0 and flag_CD19_on=0 then cell_group='CD4';
   else if flag_CD4_on=0 and flag_CD8_on=1 and flag_CD19_on=0 then cell_group='CD8';
   else if flag_CD4_on=0 and flag_CD8_on=0 and flag_CD19_on=1 then cell_group='CD19';
   else if flag_CD4_on=1 and flag_CD8_on=1 and flag_CD19_on=0 then cell_group='CD4_CD8';
   else if flag_CD4_on=1 and flag_CD8_on=0 and flag_CD19_on=1 then cell_group='CD4_CD19';
   else if flag_CD4_on=0 and flag_CD8_on=1 and flag_CD19_on=1 then cell_group='CD8_CD19';
   else if flag_CD4_on=1 and flag_CD8_on=1 and flag_CD19_on=1 then cell_group='CD4_CD8_CD19';
   else cell_group='None';
   if flag_diabetes_gene=1 then output;
   keep fusion_id gene_id flag_CD4_on flag_CD8_on flag_CD19_on cell_group;
run;

proc freq data=t1d_fusions_expressed;
   tables cell_group;
run;


data ai_fusions_expressed;
   set results_by_fusion_w_flags;
   length cell_group $12.;
   if flag_CD4_on=1 and flag_CD8_on=0 and flag_CD19_on=0 then cell_group='CD4';
   else if flag_CD4_on=0 and flag_CD8_on=1 and flag_CD19_on=0 then cell_group='CD8';
   else if flag_CD4_on=0 and flag_CD8_on=0 and flag_CD19_on=1 then cell_group='CD19';
   else if flag_CD4_on=1 and flag_CD8_on=1 and flag_CD19_on=0 then cell_group='CD4_CD8';
   else if flag_CD4_on=1 and flag_CD8_on=0 and flag_CD19_on=1 then cell_group='CD4_CD19';
   else if flag_CD4_on=0 and flag_CD8_on=1 and flag_CD19_on=1 then cell_group='CD8_CD19';
   else if flag_CD4_on=1 and flag_CD8_on=1 and flag_CD19_on=1 then cell_group='CD4_CD8_CD19';
   else cell_group='None';
   if flag_diabetes_gene=0 and flag_immuno_gene=1 then output;
   keep fusion_id gene_id flag_CD4_on flag_CD8_on flag_CD19_on cell_group;
run;

proc freq data=ai_fusions_expressed;
   tables cell_group;
run;



/* GENE COUNTS */


/* T1D GENES */
/* Split genes on cell type */

data CD4_genes;
   set results_by_fusion_w_flags;
   if flag_CD4_on=1 and flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=1;
   keep gene_id;
run;

data CD8_genes;
   set results_by_fusion_w_flags;
   if flag_CD8_on=1 and flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=1;
   keep gene_id;
run;

data CD19_genes;
   set results_by_fusion_w_flags;
   if flag_CD19_on=1 and flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=1;
   keep gene_id;
run;

data all_genes;
   set results_by_fusion_w_flags;
   if flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=1;
   keep gene_id;
run;

proc sort data=CD4_genes nodup;
  by gene_id;
run;

proc sort data=CD8_genes nodup;
  by gene_id;
run;

proc sort data=CD19_genes nodup;
  by gene_id;
run;

proc sort data=all_genes nodup;
  by gene_id;
run;

data expressed_genes_by_group_t1d;
   merge CD4_genes (in=in1) CD8_genes (in=in2) CD19_genes (in=in3) all_genes (in=in4);
   by gene_id;
   length cell_group $12.;
   if in1 and in2 and in3 then cell_group='CD4_CD8_CD19';
 
   else if in1 and not in2 and not in3 then cell_group='CD4';
   else if not in1 and in2 and not in3 then cell_group='CD8';
   else if not in1 and not in2 and in3 then cell_group='CD19';

   else if in1 and in2 and not in3 then cell_group='CD4_CD8';
   else if in1 and not in2 and in3 then cell_group='CD4_CD19';
   else if not in1 and in2 and in3 then cell_group='CD8_CD19';

   else cell_group='OFF';
run;

proc freq data=expressed_genes_by_group_t1d;
   tables cell_group;
run;


/* IMMUNO GENES */
/* Split genes on cell type */

data CD4_genes;
   set results_by_fusion_w_flags;
   if flag_CD4_on=1 and flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=0 and flag_immuno_gene=1;
   keep gene_id;
run;

data CD8_genes;
   set results_by_fusion_w_flags;
   if flag_CD8_on=1 and flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=0 and flag_immuno_gene=1;
   keep gene_id;
run;

data CD19_genes;
   set results_by_fusion_w_flags;
   if flag_CD19_on=1 and flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=0 and flag_immuno_gene=1;
   keep gene_id;
run;

data all_genes;
   set results_by_fusion_w_flags;
   if flag_multigene=0 and flag_pseudogene=0 and flag_diabetes_gene=0 and flag_immuno_gene=1;
   keep gene_id;
run;

proc sort data=CD4_genes nodup;
  by gene_id;
run;

proc sort data=CD8_genes nodup;
  by gene_id;
run;

proc sort data=CD19_genes nodup;
  by gene_id;
run;

proc sort data=all_genes nodup;
  by gene_id;
run;

data expressed_genes_by_group_ai;
   merge CD4_genes (in=in1) CD8_genes (in=in2) CD19_genes (in=in3) all_genes (in=in4);
   by gene_id;
   length cell_group $12.;
   if in1 and in2 and in3 then cell_group='CD4_CD8_CD19';
 
   else if in1 and not in2 and not in3 then cell_group='CD4';
   else if not in1 and in2 and not in3 then cell_group='CD8';
   else if not in1 and not in2 and in3 then cell_group='CD19';

   else if in1 and in2 and not in3 then cell_group='CD4_CD8';
   else if in1 and not in2 and in3 then cell_group='CD4_CD19';
   else if not in1 and in2 and in3 then cell_group='CD8_CD19';

   else cell_group='OFF';
run;

proc freq data=expressed_genes_by_group_ai;
   tables cell_group;
run;



/* FUSION COUNTS
cell_group	All	Immuno	T1D
CD4		2847	169	2
CD8		3576	223	4
CD19		7953	476	20
CD4+8		10053	795	15
CD4+19		1780	102	9
CD8+19		2384	164	2
ALL		175001	15836	283
TOTAL_ON	203594	17765	335

OFF		148365	5958	73
TOTAL_ALL	351959	23723	408
*/


/* GENE COUNTS
cell_group	All	Immuno	T1D
CD4		133	4	0
CD8		117	7	0
CD19		359	19	0
CD4+8		414	27	0
CD4+19		91	6	0
CD8+19		108	9	0
ALL		14336	1486	33
TOTAL_ON	15558	1558	33

OFF		3844	173	1
TOTAL_ALL	19402	1731	34
*/



/**** DE Fusions - each comparison ****/


/* split on comparisons */

data CD4CD8_fus_genes;
   set results_by_fusion_w_flags;
   if flag_cd4cd8_fdr05=1;
   keep fusion_id;
run;

data CD4CD19_fus_genes;
   set results_by_fusion_w_flags;
   if flag_cd4cd19_fdr05=1;
   keep fusion_id;
run;

data CD8CD19_fus_genes;
   set results_by_fusion_w_flags;
   if flag_cd8cd19_fdr05=1;
   keep fusion_id;
run;

data all_fus;
   set results_by_fusion_w_flags;
   keep fusion_id gene_id flag_pseudogene flag_immuno_gene flag_diabetes_gene;
run;


proc sort data=CD4CD8_fus_genes nodup;
   by fusion_id;
run;

proc sort data=CD4CD19_fus_genes nodup;
   by fusion_id;
run;

proc sort data=CD8CD19_fus_genes nodup;
   by fusion_id;
run;

proc sort data=all_fus nodup;
  by fusion_id;
run;

data de_fusion;
   merge CD4CD8_fus_genes (in=in1) CD4CD19_fus_genes (in=in2) CD8CD19_fus_genes (in=in3) all_fus (in=in4);
   by fusion_id;
   length comparison  $10.;
   if in1 and in2 and in3 then comparison='48_419_819';

   else if in1 and not in2 and not in3 then comparison='48';
   else if not in1 and in2 and not in3 then comparison='419';
   else if not in1 and not in2 and in3 then comparison='819';

   else if in1 and in2 and not in3 then comparison='48_419';
   else if in1 and not in2 and in3 then comparison='48_819';
   else if not in1 and in2 and in3 then comparison='419_819';

   else comparison='None';
run;


proc freq data=de_fusion;
   tables comparison;
run;


proc freq data=de_fusion;
   where flag_diabetes_gene=1;
   tables comparison;
run;

proc freq data=de_fusion;
   where flag_diabetes_gene=0 and flag_immuno_gene=1;
   tables comparison;
run;


/* GENES
COMPARISONS	TOTAL	IMMUNO	T1D
CD4/8		1331	128	2
CD4/19		10261	821	5
CD8/19		9068	728	10
CD4/8+CD4/19	6821	590	5
CD4/8+CD8/19	6452	583	18
CD4/19+CD8/19	71993	6810	79
ALL		28244	3183	107
TOTAL_ON	134170	12843	226

OFF		217789	10880	182
TOTAL		351959	23723	408
*/


/* Now do counts for genes */

data CD4CD8_DE_genes;
   set results_by_fusion_w_flags;
   if flag_cd4cd8_fdr05=1 and flag_multigene=0 and flag_pseudogene=0;
   keep gene_id;
run;

data CD4CD19_DE_genes;
   set results_by_fusion_w_flags;
   if flag_cd4cd19_fdr05=1 and flag_multigene=0 and flag_pseudogene=0;
   keep gene_id;
run;

data CD8CD19_DE_genes;
   set results_by_fusion_w_flags;
   if flag_cd8cd19_fdr05=1 and flag_multigene=0 and flag_pseudogene=0;
   keep gene_id;
run;

data all_genes;
   set results_by_fusion_w_flags;
   if flag_multigene=0 and flag_pseudogene=0;
   keep gene_id flag_pseudogene flag_immuno_gene flag_diabetes_gene;
run;


proc sort data=CD4CD8_DE_genes nodup;
   by gene_id;
run;

proc sort data=CD4CD19_DE_genes nodup;
   by gene_id;
run;

proc sort data=CD8CD19_DE_genes nodup;
   by gene_id;
run;

proc sort data=all_genes nodup;
  by gene_id;
run;

data de_genes;
   merge CD4CD8_DE_genes (in=in1) CD4CD19_DE_genes (in=in2) CD8CD19_DE_genes (in=in3) all_genes (in=in4);
   by gene_id;
   length comparison  $10.;
   if in1 and in2 and in3 then comparison='48_419_819';

   else if in1 and not in2 and not in3 then comparison='48';
   else if not in1 and in2 and not in3 then comparison='419';
   else if not in1 and not in2 and in3 then comparison='819';

   else if in1 and in2 and not in3 then comparison='48_419';
   else if in1 and not in2 and in3 then comparison='48_819';
   else if not in1 and in2 and in3 then comparison='419_819';

   else comparison='None';
run;


proc freq data=de_genes;
   tables comparison;
run;


proc freq data=de_genes;
   where flag_diabetes_gene=1;
   tables comparison;
run;

proc freq data=de_genes;
   where flag_diabetes_gene=0 and flag_immuno_gene=1;
   tables comparison;
run;


/* GENES
COMPARISONS	TOTAL	IMMUNO	T1D
CD4/8		46	5	0
CD4/19		222	17	0
CD8/19		279	233	2
CD4/8+CD4/19	293	30	0
CD4/8+CD8/19	365	33	2
CD4/19+CD8/19	5559	543	6
ALL		6976	801	22
TOTAL_ON	13740	1452	32
	
OFF		5662	279	2
TOTAL		19402	1731	34
*/



data 
   merge results_by_fusion2 (in=in1) immunogene_flags (in=in2);
   by gene_id;
   if in1;
run;


/* Enrichment of AI and T1D genes in fusion info */


data results_by_fusion_w_flags2;
   set results_by_fusion_w_flags;
   if flag_cd4cd19_fdr05=1 or flag_cd4cd8_fdr05=1 or flag_cd8cd19_fdr05=1 then flag_de=1;
   else flag_de=0;
run;

proc freq data=results_by_fusion_w_flags2;
   where flag_diabetes_gene=0;
   tables flag_de*flag_immuno_gene / chisq;
run;

/*  Table of flag_de by flag_immuno_gene

 flag_de     flag_immuno_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 183622 |  11053 | 194675
          |  59.00 |   3.55 |  62.55
          |  94.32 |   5.68 |
          |  63.95 |  45.82 |
 ---------+--------+--------+
        1 | 103496 |  13069 | 116565
          |  33.25 |   4.20 |  37.45
          |  88.79 |  11.21 |
          |  36.05 |  54.18 |
 ---------+--------+--------+
 Total      287118    24122   311240
             92.25     7.75   100.00

      Frequency Missing = 40719

            The SAS System          15:
 Statistics for Table of flag_de by flag_immuno_gene

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1   3123.1612    <.0001
Likelihood Ratio Chi-Square    1   3015.1936    <.0001
Continuity Adj. Chi-Square     1   3122.3872    <.0001
Mantel-Haenszel Chi-Square     1   3123.1511    <.0001
Phi Coefficient                       0.1002
Contingency Coefficient               0.0997
Cramer's V                            0.1002

              Fisher's Exact Test
       ----------------------------------
       Cell (1,1) Frequency (F)    183622
       Left-sided Pr <= F          1.0000
       Right-sided Pr >= F         <.0001

       Table Probability (P)       <.0001
       Two-sided Pr <= P           <.0001

         Effective Sample Size = 311240
           Frequency Missing = 40719

     WARNING: 12% of the data are missing.


*/


proc freq data=results_by_fusion_w_flags2;
   tables flag_de*flag_diabetes_gene / chisq;
run;

/* 
 Table of flag_de by flag_diabetes_gene

  flag_de     flag_diabetes_gene

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 194493 |    182 | 194675
           |  62.49 |   0.06 |  62.55
           |  99.91 |   0.09 |
           |  62.57 |  44.61 |
  ---------+--------+--------+
         1 | 116339 |    226 | 116565
           |  37.38 |   0.07 |  37.45
           |  99.81 |   0.19 |
           |  37.43 |  55.39 |
  ---------+--------+--------+
  Total      310832      408   311240
              99.87     0.13   100.00

       Frequency Missing = 40719

Statistics for Table of flag_de by flag_diabetes_gene

Statistic                     DF       Value      Prob
------------------------------------------------------
Chi-Square                     1     56.1312    <.0001
Likelihood Ratio Chi-Square    1     53.9367    <.0001
Continuity Adj. Chi-Square     1     55.3669    <.0001
Mantel-Haenszel Chi-Square     1     56.1310    <.0001
Phi Coefficient                       0.0134
Contingency Coefficient               0.0134
Cramer's V                            0.0134


 Statistics for Table of flag_de by flag_diabetes_gene

                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)    194493
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

             Effective Sample Size = 311240
               Frequency Missing = 40719

         WARNING: 12% of the data are missing.

*/


/* Enrichment on gene-level */


data de_gene_2;
  set de_genes;
  if comparison='None' then flag_de=0;
  else flag_de=1;
run;

proc freq data=de_gene_2;
   where flag_pseudogene=0 and flag_diabetes_gene=0;
   tables flag_de*flag_immuno_gene / chisq;
run;

/* 
Table of flag_de by flag_immuno_gene

flag_de     flag_immuno_gene

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   5382 |    280 |   5662
         |  27.74 |   1.44 |  29.18
         |  95.05 |   4.95 |
         |  30.51 |  15.87 |
---------+--------+--------+
       1 |  12256 |   1484 |  13740
         |  63.17 |   7.65 |  70.82
         |  89.20 |  10.80 |
         |  69.49 |  84.13 |
---------+--------+--------+
Total       17638     1764    19402
            90.91     9.09   100.00

  Statistics for Table of flag_de by flag_immuno_gene

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1    166.3253    <.0001
 Likelihood Ratio Chi-Square    1    185.1334    <.0001
 Continuity Adj. Chi-Square     1    165.6177    <.0001
 Mantel-Haenszel Chi-Square     1    166.3168    <.0001
 Phi Coefficient                       0.0926
 Contingency Coefficient               0.0922
 Cramer's V                            0.0926


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)      5382
           Left-sided Pr <= F          1.0000
           Right-sided Pr >= F         <.0001

           Table Probability (P)       <.0001
           Two-sided Pr <= P           <.0001

                  Sample Size = 19402


*/

proc freq data=de_gene_2;
   where flag_pseudogene=0;
   tables flag_de*flag_diabetes_gene / chisq;
run;


/* 
Table of flag_de by flag_diabetes_gene

 flag_de     flag_diabetes_gene

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   5660 |      2 |   5662
          |  29.17 |   0.01 |  29.18
          |  99.96 |   0.04 |
          |  29.22 |   5.88 |
 ---------+--------+--------+
        1 |  13708 |     32 |  13740
          |  70.65 |   0.16 |  70.82
          |  99.77 |   0.23 |
          |  70.78 |  94.12 |
 ---------+--------+--------+
 Total       19368       34    19402
             99.82     0.18   100.00

 Statistics for Table of flag_de by flag_diabetes_gene

 Statistic                     DF       Value      Prob
 ------------------------------------------------------
 Chi-Square                     1      8.9474    0.0028
 Likelihood Ratio Chi-Square    1     11.8134    0.0006
 Continuity Adj. Chi-Square     1      7.8536    0.0051
 Mantel-Haenszel Chi-Square     1      8.9469    0.0028
 Phi Coefficient                       0.0215
 Contingency Coefficient               0.0215
 Cramer's V                            0.0215


                  Fisher's Exact Test
           ----------------------------------
           Cell (1,1) Frequency (F)      5660
           Left-sided Pr <= F          0.9999
           Right-sided Pr >= F         0.0009

           Table Probability (P)       0.0008
           Two-sided Pr <= P           0.0019

                  Sample Size = 19402



*/
