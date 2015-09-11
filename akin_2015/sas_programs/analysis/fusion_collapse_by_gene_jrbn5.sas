
/* import libraries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';

/* Collapse fusions on gene */

proc sort data=sugrue.results_by_fusions_final;
   by fusion_id;
run;

/* drop multigene fusions */

data fusions_nomulti sugrue.multigene_fusions;
   set sugrue.results_by_fusions_final;
   if flag_multigene=0 then output fusions_nomulti;
   else output sugrue.multigene_fusions;
run;

/* get collapsed gene info */

data fusion2gene_nomulti;
    set fusions_nomulti;
    exon_count=1;
run;


data all_gene_info;
   set fusion2gene_nomulti;
   keep gene_id region_length flag_control_on flag_treat_on flag_all_on flag_p05 flag_fdr_05 flag_fdr_05 exon_count;
   run;

proc sort data=all_gene_info;
  by gene_id;
run;

proc means data=all_gene_info noprint; 
   var exon_count;
   by gene_id;
   output out=gene_info_collapsed sum=;
run;


/* Redo but for only fusions expressed in both */

proc means data=all_gene_info noprint;
   var region_length flag_control_on flag_treat_on flag_all_on flag_p05 flag_fdr_05 exon_count;
   by gene_id ;
   where flag_all_on=1;
   output out=gene_info_collapsed_exp sum=;
run;

proc means data=all_gene_info noprint;
   var region_length flag_control_on flag_treat_on flag_all_on flag_p05 flag_fdr_05 exon_count;
   by gene_id ;
   where flag_control_on=1 and flag_treat_on=0;
   output out=gene_info_collapsed_con sum=;
run;

proc means data=all_gene_info noprint;
   var region_length flag_control_on flag_treat_on flag_all_on flag_p05 flag_fdr_05 exon_count;
   by gene_id ;
   where flag_control_on=0 and flag_treat_on=1;
   output out=gene_info_collapsed_treat sum=;
run;

/* back in total exon count by gene */

data gene_info_collapsed2;
   set gene_info_collapsed;
   drop _TYPE_ _FREQ_;
   rename exon_count=total_exon_count;
run;


proc sort data=gene_info_collapsed2;
   by gene_id;
run;

proc sort data=gene_info_collapsed_exp;
   by gene_id;
run;

proc sort data=gene_info_collapsed_con;
   by gene_id;
run;

proc sort data=gene_info_collapsed_treat;
   by gene_id;
run;


data gene_info_collapsed_exp2 no_exp oops;
   merge gene_info_collapsed_exp (in=in1) gene_info_collapsed2 (in=in2);
   by gene_id;
   if in1 and in2 then output gene_info_collapsed_exp2;
   else if in2 then output no_exp;
   else output oops;
   rename exon_count=exp_exon_count;
   drop _TYPE_ _FREQ_;
run;


data gene_info_collapsed_con2 no_con oops;
   merge gene_info_collapsed_con (in=in1) gene_info_collapsed2 (in=in2);
   by gene_id;
   if in1 and in2 then output gene_info_collapsed_con2;
   else if in2 then output no_con;
   else output oops;
run;


data gene_info_collapsed_treat2 no_treat oops;
   merge gene_info_collapsed_treat (in=in1) gene_info_collapsed2 (in=in2);
   by gene_id;
   if in1 and in2 then output gene_info_collapsed_treat2;
   else if in2 then output no_treat;
   else output oops;
run;

/* total counts for exons by gene expressed in each group */

data exons_exp_both;
   set gene_info_collapsed_exp2;
   keep gene_id flag_all_on;
run;

data exons_exp_con;
   set gene_info_collapsed_con2;
   keep gene_id flag_control_on;
run;

data exons_exp_treat;
   set gene_info_collapsed_treat2;
   keep gene_id flag_treat_on;
run;

proc sort data=exons_exp_both;
   by gene_id;
run;

proc sort data=exons_exp_con;
   by gene_id;
run;

proc sort data=exons_exp_treat;
   by gene_id;
run;

data exons_by_gene_exp;
   merge exons_exp_both exons_exp_con exons_exp_treat;
   by gene_id;
   if flag_all_on=. then flag_all_on=0;
   if flag_control_on=. then flag_control_on=0;
   if flag_treat_on=. then flag_treat_on=0;
   exon_cnt_con_total=flag_all_on+flag_control_on;
   exon_cnt_treat_total=flag_all_on+flag_treat_on;
   rename flag_all_on=exon_cnt_both_only;
   rename flag_control_on=exon_cnt_con_only;
   rename flag_treat_on=exon_cnt_treat_only;
run;


/* cat up_down indicator */
/* changed this to only include differentially expressed fusions */
/* Nominal P<0.05. FDR corrected exons (0.05) only have 24 fusions */

data updown_fusions;
   set fusion2gene_nomulti;
   keep gene_id fusion_id flag_up_down;
   if flag_all_on=1 and flag_p05=1 then output;
   run;

/* get counts first */

proc freq noprint data=updown_fusions;
   tables gene_id / out=updown_cnt;
run;

proc sort data=updown_cnt;
  by descending count;
run;
*max=47 fusions per gene;

proc sort data=updown_fusions;
   by gene_id fusion_id;
run;

data updown_fusions2; 
  array updown[47] $ 1;

  retain updown1-updown47;

  set updown_fusions;
  by gene_id;
  
  if first.gene_id then do;
     call missing(of updown1-updown47);
     records = 0;
  end;

  records + 1;
  updown[records]=flag_up_down;
  if last.gene_id then output;
run;

  *clean up the output file;

data updown_fusions3;
  set updown_fusions2;
  length updown_cat $ 47;
  rename records= num_exons_cat;
         updown_cat= cats(OF updown1-updown47);
  drop updown1-updown47 flag_up_down fusion_id ;
  run;

/*merge all in and make permenant*/

proc sort data=exons_by_gene_exp;
   by gene_id;
run;

proc sort data=gene_info_collapsed_exp2;
by gene_id;
run;

proc sort data=updown_fusions3;
by gene_id;
run;


data exp_gene_info_w_updown oops;
   merge gene_info_collapsed_exp2 (in=in1) updown_fusions3 (in=in2);
   by gene_id;
   if in1 and in2 then output exp_gene_info_w_updown;
   else if in1 then do;
      num_exons_cat=0;
      updown_cat='.';
      output exp_gene_info_w_updown;
      end;
   else output oops;
   rename num_exons_cat=de_exon_cnt;
run;

data exp_gene_info_w_updown_2 no_exp oops;
   merge exp_gene_info_w_updown (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output exp_gene_info_w_updown_2;
   else if in2 then output no_exp;
   else output oops;
run;

data sugrue.results_gene_summary sugrue.non_expressed_gene oops;
   merge exp_gene_info_w_updown_2 (in=in1) gene_info_collapsed (in=in2);
   by gene_id;
   if in1 and in2 then output sugrue.results_gene_summary;
   else if in2 then output sugrue.non_expressed_gene ;
   else output oops;
   drop _TYPE_ _FREQ_;
run;

/* Make control-only and treatment-only gene summaries and make permenant */

proc sort data=gene_info_collapsed_con2;
by gene_id;
run;

proc sort data=gene_info_collapsed_treat2;
by gene_id;
run;


data sugrue.gene_info_collapsed_con no_con oops;
   merge gene_info_collapsed_con2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output sugrue.gene_info_collapsed_con;
   else if in1 then output oops;
   else output no_con;
   drop flag_p05 flag_fdr_05 _TYPE_ _FREQ_;
run;


data sugrue.gene_info_collapsed_treat no_treat oops;
   merge gene_info_collapsed_treat2 (in=in1) exons_by_gene_exp (in=in2);
   by gene_id;
   if in1 and in2 then output sugrue.gene_info_collapsed_treat;
   else if in1 then output oops;
   else output no_treat;
   drop flag_p05 flag_fdr_05 _TYPE_ _FREQ_;
run;

/* make GOI set */

data goi_exp;
   set sugrue.results_gene_summary;
      if gene_id='PAX6' then output;
      else if gene_id='FGFR2' then output;
      else if gene_id='CD44' then output;
      else if gene_id='CTNND1' then output;
      else if gene_id='TMX2andC11orf31andCTNND1' then output;
      else if gene_id='ENAH' then output;
      else if gene_id='SLC37A2' then output;
      else if gene_id='ARHGEF11' then output;
      else if gene_id='FOXJ3' then output;
      else if gene_id='FAM50A' then output;
      else if gene_id='PSENEN' then output;
      else if gene_id='PSENENandLIN37' then output;
      else if gene_id='ECT2' then output;
      else if gene_id='NCSTN' then output;
      else if gene_id='SLC1A2' then output;
      else if gene_id='NCRNA00085' then output;
      else if gene_id='NCRNA00077' then output;
      else if gene_id='HAS2AS' then output;
run;

data goi_noexp;
   set sugrue.non_expressed_gene;
      if gene_id='PAX6' then output;
      else if gene_id='FGFR2' then output;
      else if gene_id='CD44' then output;
      else if gene_id='CTNND1' then output;
      else if gene_id='TMX2andC11orf31andCTNND1' then output;
      else if gene_id='ENAH' then output;
      else if gene_id='SLC37A2' then output;
      else if gene_id='ARHGEF11' then output;
      else if gene_id='FOXJ3' then output;
      else if gene_id='FAM50A' then output;
      else if gene_id='PSENEN' then output;
      else if gene_id='PSENENandLIN37' then output;
      else if gene_id='ECT2' then output;
      else if gene_id='NCSTN' then output;
      else if gene_id='SLC1A2' then output;
      else if gene_id='NCRNA00085' then output;
      else if gene_id='NCRNA00077' then output;
      else if gene_id='HAS2AS' then output;
run;


data sugrue.gene_summary_all_goi;
   set goi_exp (in=in1) goi_noexp (in=in2);
   if in1 then flag_goi_exp=1;
   if in2 then flag_go_noexp=0;
run;


proc export data=sugrue.multigene_fusions
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/multigene_fusions_results.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.results_gene_summary
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/results_gene_summary.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.gene_summary_all_goi
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/gene_summary_all_goi.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.non_expressed_gene
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/gene_summary_noexp.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.gene_info_collapsed_con
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/gene_summary_con-only.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.gene_info_collapsed_treat
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/gene_summary_treat-only.csv'
	dbms=csv replace;
	run;

