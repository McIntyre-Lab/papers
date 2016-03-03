
/* Get counts for genes */

libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* collapse results by gene */


/* first, drop multigene fusions and pseudogenes */

data gene_results;
    set con.results_by_fusion_w_flags;
    if flag_multigene=1 then delete;
    *if flag_pseudogene=1 then delete;
    *if flag_pseudogene=. then delete;

    if flag_immuno_gene=. then flag_immuno_gene=0;
    if flag_diabetes_gene=. then flag_diabetes_gene=0;
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


proc sort data=con.immunogene_flags;
   by gene_id;
run;


data con.results_by_fusion_genes_w_flags;
   merge gene_results_all_flags2 (in=in1) con.immunogene_flags;
   by gene_id;
   if in1;
run;

/* Counts for genes */

/* All genes */
proc freq data=con.results_by_fusion_genes_w_flags noprint;
   tables flag_CD4_on*flag_CD8_on*flag_CD19_on / out=gene_counts_on_all;
run;

/*
	No pseudo	Pseudo
CD4	150		607
CD8	142		630
CD19	405		1720
CD4/8	492		1948
CD4/19	104		444	
CD8/19	120		462
ALL	15694		33070
OFF	4781		25355
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
	No pseudo	pseudo
4v8	48		155
4v19	341		1475
8v19	358		1479
48/419	352		717
48/819	418		891
419/819	6379		13011
ALL	6770		8579
NOT DE	1014		6737
NO TEST	6208		31192

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
NO TEST	247
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
NO TEST	1
*/

