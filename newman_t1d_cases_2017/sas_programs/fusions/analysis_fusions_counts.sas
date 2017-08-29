/* Get counts for fusions and genes */

libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* Merge in gene_ids with fusions */


data fusion2gene;
   set fus.unique_info_fusions_si;
   keep fusion_id gene_id flag_multigene;
run;

data fusion_results;
   set con.results_by_fusion_w_fdr;
   keep fusion_id flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
run;

data fusion_flags;
   set con.fusions_on_gt_apn0;
   keep fusion_id flag_cd19_on flag_cd8_on flag_cd4_on flag_fusion_all_on0;
run;

proc sort data=fusion2gene;
   by fusion_id;
proc sort data=fusion_results;
   by fusion_id;
proc sort data=fusion_flags;
   by fusion_id;
run;


data fusion2gene_w_flags;
   merge fusion2gene fusion_results fusion_flags;
   by fusion_id;
run;


/* Merge in Immunoflags */

data immunogene_flags;
   set con.immunogene_flags;
run;

proc sort data=immunogene_flags;
   by gene_id;
proc sort data=fusion2gene_w_flags;
   by gene_id;
run;

data fusion2gene_for_counts;
   merge immunogene_flags (in=in1) fusion2gene_w_flags (in=in2);
   by gene_id;
   if in1 and in2 then output;
   else if in1 then delete;
   else output;
run;



/* Make permenant */

data con.results_by_fusion_w_flags;
   set fusion2gene_for_counts;
run;


/* Counts for fusions */

/* All genes */
proc freq data=con.results_by_fusion_w_flags noprint;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=fusion_counts_on_all;
run;

/*
CD4	2709
CD8	3483
CD19	7756
CD4/8	9698
CD4/19	1709
CD8/19	2321
ALL	163713
OFF	138168
*/

/* Autoimmune genes */

proc freq data=con.results_by_fusion_w_flags noprint;
   where flag_immuno_gene=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=fusion_counts_on_ai;
run;

/*
CD4	167
CD8	219
CD19	477
CD4/8	781
CD4/19	107
CD8/19	155
ALL	14656
OFF	5005
*/


/* T1D genes */
proc freq data=con.results_by_fusion_w_flags noprint;
   where flag_diabetes_gene=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=fusion_counts_on_t1d;
run;

/*
CD4	1
CD8	4
CD19	20
CD4/8	15
CD4/19	8
CD8/19	2
ALL	262
OFF	62
*/


/* Counts for DE fusions */

/* All genes */
proc freq data=con.results_by_fusion_w_flags noprint;
   where flag_fusion_all_on0=1;
   tables flag_CD4CD8_fdr05*flag_CD4CD19_fdr05*flag_CD8CD19_fdr05 / out=fusion_de_all;
run;


/*
4v8	1257
4v19	9641
8v19	8488
48/419	6244
48/819	6036
419/819	67292
ALL	26201
NOT DE	38554
*/

/* Autoimmune genes */

proc freq data=con.results_by_fusion_w_flags noprint;
   where flag_fusion_all_on0=1 and flag_immuno_gene=1;
   tables flag_CD4CD8_fdr05*flag_CD4CD19_fdr05*flag_CD8CD19_fdr05 / out=fusion_de_ai;
run;


/*
4v8	134
4v19	761
8v19	678
48/419	538
48/819	557
419/819	6257
ALL	2960
NOT DE	2771

*/

/* T1D genes */
proc freq data=con.results_by_fusion_w_flags noprint;
   where flag_fusion_all_on0=1 and flag_diabetes_gene=1;
   tables flag_CD4CD8_fdr05*flag_CD4CD19_fdr05*flag_CD8CD19_fdr05 / out=fusion_de_t1d;
run;



/*
4v8	3
4v19	4
8v19	8
48/419	5
48/819	18
419/819	74
ALL	98
NOT DE	52
*/
