/* import libraries */

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data';


/* Get covariates: sex pool V3 subject_id */

data covar_for_anova_&sysparm.;
   set mysas.all_covariates;
   keep Name cell_type subject_id pool sex V3;
run;


/* subset expressed genes for ANOVAs */

data junc_info_&sysparm.;
   set mysas.counts_by_splicing_w_flags_&sysparm.;
   keep event_id flag_CD19_on flag_CD4_on flag_CD8_on flag_all_on;
run;

proc sort data=junc_info_&sysparm. nodup;
   by event_id;
run;

proc sort data=mysas.splicing_counts_w_zeros_&sysparm.;
   by event_id;

data splicing_counts_w_flags_&sysparm.;
   merge junc_info_&sysparm. (in=in1) mysas.splicing_counts_w_zeros_&sysparm. (in=in2);
   by event_id;
   if in1 and in2 then output;
run;


data splicing_all_on_&sysparm.;
   set splicing_counts_w_flags_&sysparm.;
   log_depth=log(depth + 1);
   if flag_all_on=1 then output;
run;

/* merge in covariates */

proc sort data=splicing_all_on_&sysparm.;
    by Name;
run;

proc sort data=covar_for_anova_&sysparm.;
   by Name;
run;

data splicing_all_on_w_cov_&sysparm. oops low_cov_sample;
   merge covar_for_anova_&sysparm. (in=in1) splicing_all_on_&sysparm. (in=in2);
   by Name;
   if in1 and in2 then output splicing_all_on_w_cov_&sysparm.;
   else if in1 then output oops;
   else output low_cov_sample;
run;


/* ANOVAs */

proc sort data=splicing_all_on_w_cov_&sysparm.;
  by event_id Name;
  run;


ods listing close;
proc glimmix data=splicing_all_on_w_cov_&sysparm.;
   by event_id;
   class cell_type sex pool subject_id;
   model log_depth = cell_type | sex pool v3 / ddfm=kr;
   random pool;
   random resid / subject=subject_id;

output out=splicing_resid_&sysparm. resid=resid pred=pred student=stu;


lsmeans cell_type/diff;			/*CD19		CD4		CD8*/
contrast 'CD4-CD8' cell_type		0		 1		-1	;
contrast 'CD4-CD19' cell_type		1		-1	 	 0	;
contrast 'CD8-CD19' cell_type	 	1		 0		-1	;

estimate 'CD4-CD8' cell_type		0		 1		-1	;
estimate 'CD4-CD19' cell_type		1		-1		 0	;
estimate 'CD8-CD19' cell_type		1		 0		-1	;

ods output tests3=splicing_anova_&sysparm.
	lsmeans = splicing_lsmeans_&sysparm.
	diffs = splicing_diffs_&sysparm.
	contrasts = splicing_contrasts_&sysparm.
	estimates = splicing_estimates_&sysparm.;


run;
quit;

* Flag residuals;
proc univariate data = splicing_resid_&sysparm. normal noprint;
  by event_id;
  var Resid;
  output out = splicing_normtest_&sysparm. probn=pnorm;
  run;
data splicing_flag_resids_&sysparm.;
  set splicing_normtest_&sysparm.;
  if pnorm = . then flag_fail_norm = .;
  	else if pnorm le 0.05 then flag_fail_norm = 1;
	else flag_fail_norm = 0;
  run;


/* Split contrasts and remerge */

data cd4cd8_contrast_&sysparm. cd4cd19_contrast_&sysparm. cd8cd19_contrast_&sysparm.;
   set splicing_contrasts_&sysparm.;
   if label='CD4-CD8' then output cd4cd8_contrast_&sysparm.;
   if label='CD4-CD19' then output cd4cd19_contrast_&sysparm.;
   if label='CD8-CD19' then output cd8cd19_contrast_&sysparm.;
   keep event_id ProbF;
run;

data cd4cd8_contrast2_&sysparm.;
    set cd4cd8_contrast_&sysparm.;
    if ProbF lt 0.05 then flag_p05_cd4cd8=1;
    else flag_p05_cd4cd8=0;
    rename ProbF=CD4_CD8_P;
run;

data cd4cd19_contrast2_&sysparm.;
    set cd4cd19_contrast_&sysparm.;
    if ProbF lt 0.05 then flag_p05_cd4cd19=1;
    else flag_p05_cd4cd19=0;
    rename ProbF=CD4_CD19_P;
run;

data cd8cd19_contrast2_&sysparm.;
    set cd8cd19_contrast_&sysparm.;
    if ProbF lt 0.05 then flag_p05_cd8cd19=1;
    else flag_p05_cd8cd19=0;
    rename ProbF=CD8_CD19_P;
run;

proc sort data=cd4cd8_contrast2_&sysparm.;
   by event_id;
run;

proc sort data=cd4cd19_contrast2_&sysparm.;
   by event_id;
run;

proc sort data=cd8cd19_contrast2_&sysparm.;
   by event_id;
run;

data contrasts_merge_&sysparm.;
   merge cd4cd8_contrast2_&sysparm. cd4cd19_contrast2_&sysparm. cd8cd19_contrast2_&sysparm.;
   by event_id;
run;

data splicing_anova_2_&sysparm.;
   set splicing_anova_&sysparm.;
   if effect='cell_type';
run;




/* merge results */

proc sort data=contrasts_merge_&sysparm.;
by event_id;
run;

proc sort data=splicing_flag_resids_&sysparm.;
by event_id;
run;

proc sort data=splicing_anova_2_&sysparm.;
    by event_id;
run;



data mysas.splicing_anova_results_&sysparm.;
   merge contrasts_merge_&sysparm. (in=in1) splicing_flag_resids_&sysparm. (in=in2) splicing_anova_2_&sysparm. (in=in3);
   by event_id;
   if in1 and in2 and in3 then output;
run;

