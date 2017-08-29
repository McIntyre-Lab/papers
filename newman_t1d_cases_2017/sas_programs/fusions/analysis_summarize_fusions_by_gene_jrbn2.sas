/*** IMPORT FUSION DATA ***/

/* import libraries */
libname con '/home/jrbnewman/concannon/sas_data';


/* Collapse fusions on gene */

proc sort data=con.results_by_fusion_final;
   by fusion_id;
run;

/* drop multigene fusions */

data fusions_nomulti con.multigene_fusions;
   set con.results_by_fusion_final;
   if flag_multigene=0 then output fusions_nomulti;
   else output con.multigene_fusions;
run;

/* get collapsed gene info */

data fusion2gene_nomulti;
    set fusions_nomulti;
    exon_count=1;
run;

data all_gene_info;
   set fusion2gene_nomulti;
   keep gene_id flag_CD19_on flag_CD4_on flag_CD8_on  flag_fusion_all_on0 flag_p05_cd4_cd8 flag_cd4cd8_fdr05 flag_p05_cd4_cd19 flag_cd4cd19_fdr05 flag_p05_cd8_cd19 flag_cd8cd19_fdr05 exon_count;
   run;

proc sort data=all_gene_info;
  by gene_id;
run;

/* Calculate total number of exons per gene */
proc means data=all_gene_info noprint; 
   var exon_count;
   by gene_id;
   output out=gene_info_collapsed sum=;
run;


/* Calculate number of fusions expressed per gene */

* fusions all on;
proc means data=all_gene_info noprint;
   var flag_fusion_all_on0 flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
   by gene_id ;
   where flag_fusion_all_on0=1;
   output out=gene_info_collapsed_exp sum=;
run;

* CD19 only;
proc means data=all_gene_info noprint;
   var exon_count;
   by gene_id ;
   where flag_CD19_on=1 and flag_CD4_on=0 and flag_CD8_on=0;
   output out=gene_info_collapsed_CD19 sum=flag_cd19_on;
run;

* CD4 only;
proc means data=all_gene_info noprint;
   var exon_count;
   by gene_id ;
   where flag_CD19_on=0 and flag_CD4_on=1 and flag_CD8_on=0;
   output out=gene_info_collapsed_CD4 sum=flag_cd4_on;
run;

* CD8 only;
proc means data=all_gene_info noprint;
   var exon_count;
   by gene_id ;
   where flag_CD19_on=0 and flag_CD4_on=0 and flag_CD8_on=1;
   output out=gene_info_collapsed_CD8 sum=flag_cd8_on;
run;

* fusions on in only CD4 and CD19;
proc means data=all_gene_info noprint;
   var exon_count;
   by gene_id ;
   where flag_CD19_on=1 and flag_CD4_on=1 and flag_CD8_on=0;
   output out=gene_info_collapsed_CD4CD19 sum=flag_cd4cd19_on;
run;

* fusions on in only CD8 and CD19;
proc means data=all_gene_info noprint;
   var exon_count;
   by gene_id ;
   where flag_CD19_on=1 and flag_CD4_on=0 and flag_CD8_on=1;
   output out=gene_info_collapsed_CD8CD19 sum=flag_cd8cd19_on;
run;

* fusions on in only CD4 and CD8;
proc means data=all_gene_info noprint;
   var exon_count;
   by gene_id ;
   where flag_CD19_on=0 and flag_CD4_on=1 and flag_CD8_on=1;
   output out=gene_info_collapsed_CD4CD8 sum=flag_cd4cd8_on;
run;


/* total counts for exons by gene expressed in each group */

data gene_info_collapsed2;
   set gene_info_collapsed;
   keep gene_id exon_count;
run;

data exons_exp_both;
   set gene_info_collapsed_exp;
   keep gene_id flag_fusion_all_on0 flag_cd4cd8_fdr05 flag_cd4cd19_fdr05 flag_cd8cd19_fdr05;
run;

data exons_exp_CD19;
   set gene_info_collapsed_CD19;
   keep gene_id flag_CD19_on;
run;

data exons_exp_CD4;
   set gene_info_collapsed_CD4;
   keep gene_id flag_CD4_on;
run;

data exons_exp_CD8;
   set gene_info_collapsed_CD8;
   keep gene_id flag_CD8_on;
run;

data exons_exp_CD4CD8;
   set gene_info_collapsed_CD4CD8;
   keep gene_id flag_cd4cd8_on;
run;

data exons_exp_CD4CD19;
   set gene_info_collapsed_CD4CD19;
   keep gene_id flag_CD4CD19_on;
run;

data exons_exp_CD8CD19;
   set gene_info_collapsed_CD8CD19;
   keep gene_id flag_CD8CD19_on;
run;


proc sort data=exons_exp_both;
   by gene_id;
run;

proc sort data=exons_exp_CD19;
   by gene_id;
run;

proc sort data=exons_exp_CD4;
   by gene_id;
run;

proc sort data=exons_exp_CD8;
   by gene_id;
run;

proc sort data=exons_exp_CD4CD8;
   by gene_id;
run;

proc sort data=exons_exp_CD4CD19;
   by gene_id;
run;

proc sort data=exons_exp_CD8CD19;
   by gene_id;
run;




data exons_by_gene_exp;
   merge exons_exp_both exons_exp_CD19 exons_exp_CD4 exons_exp_CD8 exons_exp_CD4CD8 exons_exp_CD4CD19 exons_exp_CD8CD19;
   by gene_id;
   if flag_fusion_all_on0=. then flag_fusion_all_on0=0;
   if flag_CD19_on=. then flag_CD19_on=0;
   if flag_CD4_on=. then flag_CD4_on=0;
   if flag_CD8_on=. then flag_CD8_on=0;
   if flag_CD4CD8_on=. then flag_CD4CD8_on=0;
   if flag_CD4CD19_on=. then flag_CD4CD19_on=0;
   if flag_CD8CD19_on=. then flag_CD8CD19_on=0;

   if flag_CD4CD8_fdr05=. then flag_CD4CD8_fdr05=0;
   if flag_CD4CD19_fdr05=. then flag_CD4CD19_fdr05=0;
   if flag_CD8CD19_fdr05=. then flag_CD8CD19_fdr05=0;

   exon_cnt_CD19_total=flag_fusion_all_on0+flag_CD19_on+flag_CD4CD19_on+flag_CD8CD19_on;
   exon_cnt_CD4_total=flag_fusion_all_on0+flag_CD4_on+flag_CD4CD19_on+flag_CD4CD8_on;
   exon_cnt_CD8_total=flag_fusion_all_on0+flag_CD8_on+flag_CD4CD8_on+flag_CD8CD19_on;

   rename flag_fusion_all_on0=exon_cnt_all_only;
   rename flag_CD19_on=exon_cnt_CD19_only;
   rename flag_CD4_on=exon_cnt_CD4_only;
   rename flag_CD8_on=exon_cnt_CD8_only;
   rename flag_CD4CD8_on=exon_cnt_CD4_CD8_only;
   rename flag_CD4CD19_on=exon_cnt_CD4_CD19_only;
   rename flag_CD8CD19_on=exon_cnt_CD8_CD19_only;

run;


/*merge all in and make permenant*/

proc sort data=exons_by_gene_exp;
   by gene_id;
run;

proc sort data=gene_info_collapsed2;
by gene_id;
run;


data exp_gene_info_w_exons no_exp;
   merge gene_info_collapsed2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output exp_gene_info_w_exons;
   else output no_exp;
run;


data con.results_gene_summary2 con.non_expressed_gene2 oops;
   merge exp_gene_info_w_exons (in=in1) gene_info_collapsed (in=in2);
   by gene_id;
   if in1 and in2 then output con.results_gene_summary2 ;
   else if in2 then output con.non_expressed_gene2 ;
   else output oops;
   drop _TYPE_ _FREQ_;
run;

/* Make cell type specific gene summaries and make permenant */

proc sort data=gene_info_collapsed_CD19_2;
by gene_id;
run;

proc sort data=gene_info_collapsed_CD4_2;
by gene_id;
run;


proc sort data=gene_info_collapsed_CD8_2;
by gene_id;
run;


proc sort data=gene_info_collapsed_CD4CD19_2;
by gene_id;
run;

proc sort data=gene_info_collapsed_CD4CD8_2;
by gene_id;
run;


proc sort data=gene_info_collapsed_CD8CD19_2;
by gene_id;
run;



data con.gene_info_collapsed_CD19 no_CD19 oops;
   merge gene_info_collapsed_CD19_2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output con.gene_info_collapsed_CD19;
   else if in1 then output oops;
   else output no_CD19;
   drop flag_p05_CD4_CD8 flag_cd4cd8_fdr05 flag_p05_CD4_CD19 flag_cd4cd19_fdr05 flag_p05_CD8_CD19 flag_cd8cd19_fdr05 _TYPE_ _FREQ_;
run;

data con.gene_info_collapsed_CD4 no_CD4 oops;
   merge gene_info_collapsed_CD4_2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output con.gene_info_collapsed_CD4;
   else if in1 then output oops;
   else output no_CD4;
   drop flag_p05_CD4_CD8 flag_cd4cd8_fdr05 flag_p05_CD4_CD19 flag_cd4cd19_fdr05 flag_p05_CD8_CD19 flag_cd8cd19_fdr05 _TYPE_ _FREQ_;
run;

data con.gene_info_collapsed_CD8 no_CD8 oops;
   merge gene_info_collapsed_CD8_2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output con.gene_info_collapsed_CD8;
   else if in1 then output oops;
   else output no_CD8;
   drop flag_p05_CD4_CD8 flag_cd4cd8_fdr05 flag_p05_CD4_CD19 flag_cd4cd19_fdr05 flag_p05_CD8_CD19 flag_cd8cd19_fdr05 _TYPE_ _FREQ_;
run;


data con.gene_info_collapsed_CD4CD8 no_CD4CD8 oops;
   merge gene_info_collapsed_CD4CD8_2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output con.gene_info_collapsed_CD4CD8;
   else if in1 then output oops;
   else output no_CD4CD8;
   drop flag_p05_CD4_CD8 flag_cd4cd8_fdr05 flag_p05_CD4_CD19 flag_cd4cd19_fdr05 flag_p05_CD8_CD19 flag_cd8cd19_fdr05 _TYPE_ _FREQ_;
run;



data con.gene_info_collapsed_CD4CD19 no_CD4CD19 oops;
   merge gene_info_collapsed_CD4CD19_2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output con.gene_info_collapsed_CD4CD19;
   else if in1 then output oops;
   else output no_CD4CD19;
   drop flag_p05_CD4_CD8 flag_cd4cd8_fdr05 flag_p05_CD4_CD19 flag_cd4cd19_fdr05 flag_p05_CD8_CD19 flag_cd8cd19_fdr05 _TYPE_ _FREQ_;
run;


data con.gene_info_collapsed_CD8CD19 no_CD8CD19 oops;
   merge gene_info_collapsed_CD8CD19_2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output con.gene_info_collapsed_CD8CD19;
   else if in1 then output oops;
   else output no_CD8CD19;
   drop flag_p05_CD4_CD8 flag_cd4cd8_fdr05 flag_p05_CD4_CD19 flag_cd4cd19_fdr05 flag_p05_CD8_CD19 flag_cd8cd19_fdr05 _TYPE_ _FREQ_;
run;





proc export data=con.multigene_fusions
	outfile='/home/jrbnewman/concannon/generated_files/multigene_fusions_results.csv'
	dbms=csv replace;
	run;

proc export data=con.results_gene_summary
	outfile='/home/jrbnewman/concannon/generated_files/results_gene_summary.csv'
	dbms=csv replace;
	run;

proc export data=con.non_expressed_gene
	outfile='/home/jrbnewman/concannon/generated_files/gene_summary_noexp.csv'
	dbms=csv replace;
	run;

proc export data=con.gene_info_collapsed_CD19
	outfile='/home/jrbnewman/concannon/generated_files/gene_summary_CD19-only.csv'
	dbms=csv replace;
	run;

proc export data=con.gene_info_collapsed_CD4
	outfile='/home/jrbnewman/concannon/generated_files/gene_summary_CD4-only.csv'
	dbms=csv replace;
	run;

proc export data=con.gene_info_collapsed_CD8
	outfile='/home/jrbnewman/concannon/generated_files/gene_summary_CD8-only.csv'
	dbms=csv replace;
	run;

proc export data=con.gene_info_collapsed_CD4CD8
	outfile='/home/jrbnewman/concannon/generated_files/gene_summary_CD4_CD8.csv'
	dbms=csv replace;
	run;

proc export data=con.gene_info_collapsed_CD4CD19
	outfile='/home/jrbnewman/concannon/generated_files/gene_summary_CD4_CD19.csv'
	dbms=csv replace;
	run;

proc export data=con.gene_info_collapsed_CD8CD19
	outfile='/home/jrbnewman/concannon/generated_files/gene_summary_CD8_CD19.csv'
	dbms=csv replace;
	run;


