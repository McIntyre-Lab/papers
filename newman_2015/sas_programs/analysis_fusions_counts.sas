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



/* Now collapse results by gene */


/* first, drop multigene fusions and pseudogenes */

data gene_results;
    set con.results_by_fusion_w_flags;
    if flag_multigene=1 then delete;
    if flag_pseudogene=1 then delete;
    if flag_pseudogene=. then delete;
run;


/* Sum on and de flags */

proc sort data=gene_results;
   by gene_id;
run;

proc means data=gene_results noprint;
   by gene_id;
   var flag_cd4_on;
   output out=gene_results_cd4_on sum=;
run;

proc means data=gene_results noprint;
   by gene_id;
   var flag_cd8_on;
   output out=gene_results_cd8_on sum=;
run;

proc means data=gene_results noprint;
   by gene_id;
   var flag_cd19_on;
   output out=gene_results_cd19_on sum=;
run;


proc means data=gene_results noprint;
   where flag_fusion_all_on0=1;
   by gene_id;
   var flag_cd4cd8_fdr05;
   output out=gene_results_cd48_de sum=;
run;


proc means data=gene_results noprint;
   where flag_fusion_all_on0=1;
   by gene_id;
   var flag_cd4cd19_fdr05;
   output out=gene_results_cd419_de sum=;
run;


proc means data=gene_results noprint;
   where flag_fusion_all_on0=1;
   by gene_id;
   var flag_cd8cd19_fdr05;
   output out=gene_results_cd819_de sum=;
run;


proc sort data=gene_results_cd4_on;
   by gene_id;
proc sort data=gene_results_cd8_on;
   by gene_id;
proc sort data=gene_results_cd19_on;
   by gene_id;
proc sort data=gene_results_cd48_de;
   by gene_id;
proc sort data=gene_results_cd419_de;
   by gene_id;
proc sort data=gene_results_cd819_de;
   by gene_id;
run;


data gene_results_all_flags;
   merge gene_results_cd4_on gene_results_cd8_on gene_results_cd19_on gene_results_cd48_de gene_results_cd419_de gene_results_cd819_de;
   by gene_id;
   drop _TYPE_ _FREQ_;
run;


data gene_results_all_flags2;
   set gene_results_all_flags;
   if flag_cd4_on ge 1 then flag_cd4_on=1;
   if flag_cd8_on ge 1 then flag_cd8_on=1;
   if flag_cd19_on ge 1 then flag_cd19_on=1;
   if flag_cd4cd8_fdr05 ge 1 then flag_cd4cd8_fdr05=1;
   if flag_cd4cd19_fdr05 ge 1 then flag_cd4cd19_fdr05=1;
   if flag_cd8cd19_fdr05 ge 1 then flag_cd8cd19_fdr05=1;
run;


proc sort data=gene_results_all_flags2;
   by gene_id;
run;

proc sort data=immunogene_flags;
   by gene_id;
run;


data con.results_by_fusion_genes_w_flags;
   merge gene_results_all_flags2 (in=in1) immunogene_flags;
   by gene_id;
   if in1;
run;

/* Counts for genes */

/* All genes */
proc freq data=con.results_by_fusion_genes_w_flags noprint;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=gene_counts_on_all;
run;

/*
CD4	150
CD8	142
CD19	405
CD4/8	492
CD4/19	104
CD8/19	120
ALL	15694
OFF	4781
*/

/* Autoimmune genes */
proc freq data=con.results_by_fusion_genes_w_flags noprint;
   where flag_immuno_gene=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=gene_counts_on_ai;
run;

/*
CD4	7
CD8	6
CD19	17
CD4/8	29
CD4/19	5
CD8/19	8
ALL	1488
OFF	174
*/

/* T1D genes */
proc freq data=con.results_by_fusion_genes_w_flags noprint;
   where flag_diabetes_gene=1;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=gene_counts_on_t1d;
run;

/*
CD4	0
CD8	0
CD19	0
CD4/8	0
CD4/19	0
CD8/19	0
ALL	33
OFF	1
*/




/* All genes */
proc freq data=con.results_by_fusion_genes_w_flags noprint;
   tables flag_CD4CD8_fdr05*flag_CD4CD19_fdr05*flag_CD8CD19_fdr05 / out=gene_de_all;
run;


/*
4v8	48
4v19	341
8v19	1014
48/419	352
48/819	418
419/819	6379
ALL	6770
NOT DE	1014

*/

/* Autoimmune genes */
proc freq data=con.results_by_fusion_genes_w_flags noprint;
   where flag_immuno_gene=1;
   tables flag_CD4CD8_fdr05*flag_CD4CD19_fdr05*flag_CD8CD19_fdr05 / out=gene_de_ai;
run;


/*
4v8	5
4v19	15
8v19	27
48/419	31
48/819	35
419/819	545
ALL	788
NOT DE	41
*/


/* T1D genes */
proc freq data=con.results_by_fusion_genes_w_flags noprint;
   where flag_diabetes_gene=1;
   tables flag_CD4CD8_fdr05*flag_CD4CD19_fdr05*flag_CD8CD19_fdr05 / out=gene_de_t1d;
run;


/*
4v8	0
4v19	0
8v19	2
48/419	0
48/819	2
419/819	6
ALL	22
NOT DE	1
*/

