/* import libraries */
libname fusion '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname sugrue '/home/jrbnewman/McLab/sugrue/sas_data';


/* Merge in all! */

/* counts_w_means, fusion_data_flip, anova_logapn (HypothesisType=3), results_by_fusion_w_fdr2 */

/* Means */
data counts_w_means;
   set sugrue.counts_w_means;
   keep fusion_id flag_control_on flag_treat_on flag_all_on mean_apn_control mean_apn_treat
   flag_up_down fold_change flag_low_exp_con flag_low_exp_treat flag_low_exp_both;
run;

proc sort data=counts_w_means nodup;
   by fusion_id;
run;

/* Counts */
data fusion_counts_by_subject;
   set sugrue.fusion_data_flip;
   drop _NAME_;
run;

proc sort data=fusion_counts_by_subject nodup;
   by fusion_id;
run;

/* ANOVA results */

data anova_results;
    set sugrue.anova_logapn;
    where HypothesisType=3;
    keep fusion_id FValue ProbF;
run;

proc sort data=anova_results;
    by fusion_id;
run;


data fdr_results;
  set sugrue.results_by_fusion_w_fdr2;
  keep fusion_id fdr_p flag_fdr_05;
run;

proc sort data=fdr_results;
    by fusion_id;
run;


data counts_and_means oops1 oops2;
   merge counts_w_means (in=in1) fusion_counts_by_subject (in=in2);
   by fusion_id;
   if in1 and in2 then output counts_and_means;
   else if in1 then output oops1;
   else output oops2;
run;



data counts_means_anova oops;
   merge counts_and_means (in=in1) anova_results (in=in2);
   by fusion_id;
   if in1 and in2 then output counts_means_anova;
   else if in1 then output counts_means_anova;
   else output oops;
run;


data counts_means_anova_fdr oops;
   merge counts_means_anova (in=in1) fdr_results (in=in2);
   by fusion_id;
   if in1 and in2 then output counts_means_anova_fdr;
   else if in1 then output counts_means_anova_fdr;
   else output oops;
run;



/* add in fusion2gene info */

proc sort data=fusion.unique_info_fusions_sd;
    by fusion_id;
run;

proc sort data=sugrue.results_by_fusion_w_fdr2;
    by fusion_id;
run;

data results_by_fusion_w_info oops1 oops2;
   merge fusion.unique_info_fusions_sd (in=in1) sugrue.results_by_fusion_w_fdr2 (in=in2);
   by fusion_id;
   if in1 and in2 then output results_by_fusion_w_info;
   else if in1 then output oops1;
   else output oops2;
run;


/* add in flag_sugrue_goi */

/* 
Pax6				PAX6
FGFr2				FGFR2
CD44				CD44
CtnnD1				CTNND1 or TMX2andC11orf31andCTNND1
Enah				ENAH
Slc37A2				SLC37A2
ARHGEF11			ARHGEF11
FoxJ3				FOXJ3
Fam50A				FAM50A
Psenen (Pen2)			PSENEN or PSENENandLIN37
ECT2				ECT2
NCSTN				NCSTN
GLT				SLC1A2
Linc00085/ SPACA6P		NCRNA00085
HAS2-AS1			NCRNA00077 or HAS2AS

*/

data results_fusions_goi;
   set results_by_fusion_w_info;
      if gene_id='PAX6'
      or gene_id='FGFR2'
      or gene_id='CD44'
      or gene_id='CTNND1'
      or gene_id='TMX2andC11orf31andCTNND1'
      or gene_id='ENAH'
      or gene_id='ARHGEF11'
      or gene_id='FOXJ3'
      or gene_id='FAM50A'
      or gene_id='PSENEN'
      or gene_id='PSENENandLIN37'
      or gene_id='ECT2'
      or gene_id='NCSTN'
      or gene_id='SLC1A2'
      or gene_id='NCRNA00085'
      or gene_id='NCRNA00077'
      or gene_id='HAS2AS' then flag_goi=1;
      else flag_goi=0;
run;

/* Make permenant */


data sugrue.results_by_fusions_final;
set results_fusions_goi;
drop dependent hypothesistype source ss ms ;
run;



/* split fusions on expression */

data sugrue.counts_all_goi sugrue.counts_all_on counts_con_only counts_treat_only counts_all_off oops;
   set sugrue.results_by_fusions_final;
   if flag_goi=1 then output sugrue.counts_all_goi;
   if flag_control_on=1 and flag_treat_on=1 then output sugrue.counts_all_on;
   else if flag_control_on=1 and flag_treat_on=0 then output counts_con_only;
   else if flag_control_on=0 and flag_treat_on=1 then output counts_treat_only;
   else if flag_control_on=0 and flag_treat_on=0 then output counts_all_off;
   else output oops;
run;

/* Make control-only, treatment-only and off-fusions permenant */

data sugrue.counts_con_only;
   set counts_con_only;
   drop flag_control_on flag_treat_on flag_all_on DF Fvalue ProbF pnorm flag_fail_norm flag_p05 fdr_p flag_fdr_05 flag_goi;
run;

data sugrue.counts_treat_only;
   set counts_treat_only;
   drop flag_control_on flag_treat_on flag_all_on DF Fvalue ProbF pnorm flag_fail_norm flag_p05 fdr_p flag_fdr_05 flag_goi;
run;

data sugrue.counts_all_off;
   set counts_all_off;
   drop flag_control_on flag_treat_on flag_all_on DF Fvalue ProbF pnorm flag_fail_norm flag_p05 fdr_p flag_fdr_05 flag_goi;
run;

/* export data for now */

proc export data=sugrue.counts_all_goi
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/results_fusions_goi.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.counts_all_on
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/results_fusions_all_on.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.counts_con_only
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/results_fusions_con_on_only.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.counts_treat_only
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/results_fusions_treat_on_only.csv'
	dbms=csv replace;
	run;

proc export data=sugrue.counts_all_off
	outfile='/home/jrbnewman/McLab/sugrue/pipeline_output/results_fusions_all_off.csv'
	dbms=csv replace;
	run;


