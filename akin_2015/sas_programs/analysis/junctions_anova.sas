/* JUNCTIONS */

/* Because the junctions datasets are large I will be working on this locally. I will save the final datasets to the share */

libname mysas '/scratch/lfs/sugrue/sas_analysis/sas_data';

/* subset expressed genes for ANOVAs */

proc sort data=mysas.jnc_counts_w_flags_&sysparm.;
   by fusion_id;
run;

data jnc_exp_only;
   merge mysas.jnc_counts_w_flags_&sysparm. (in=in1) mysas.junc_counts_merge_&sysparm. (in=in2);
   by fusion_id;
   if in1 and in2 then output;
run;

/* ANOVAs */

proc sort data=jnc_exp_only  nodup;
  by fusion_id subject;
  run;

/*only for checking residuals */

	proc sort data=jnc_exp_only;
	by fusion_id;
	run;

	/* freq for order */
	proc freq data=jnc_exp_only;
	tables treatment;
	run;

ods listing close;
	ods pdf close;
	ods html close;
	ods graphics off;
proc glm data=jnc_exp_only;
by fusion_id;
class treatment;
model log_apn=treatment;  

lsmeans treatment;			/*0	1*/
contrast '0-1' treatment		-1	1;
estimate '0-1' treatment		-1	1;

ods output modelanova = anova
	lsmeans = lsmean
	contrasts = contrasts
	estimates = estimates;
	output out=jnc_resid predicted=pred residual=resid student=student;
run;
quit;
ods listing;


/* Flag residuals that fail normality */
ods graphics off;
proc univariate data = jnc_resid normal noprint;
  by fusion_id;
  var Resid;
  output out = normtest probn=pnorm;
  run;


data flag_resids_logapn;
  set normtest;
  if pnorm=. then flag_fail_norm=.;
  else if pnorm ge 0 and pnorm le 0.05 then flag_fail_norm = 1;
  else flag_fail_norm = 0;
  run;

ods graphics on;
ods listing;
	ods pdf;
	ods html;

proc freq data = flag_resids_logapn noprint;
  tables flag_fail_norm / out=resid_fail_norm ;
  run;


/* flag fail norm summary 

29192 junctions passed normality (93.53%)
2018 junctions failed normality (6.47%)

*/ 


data mysas.junc_anova_&sysparm. ;
set anova ;
run ;
data mysas.junc_lsmean_&sysparm. ;
set lsmean ;
run ;

data mysas.junc_resid_flags_&sysparm.;
set flag_resids_logapn;
run;



/* quick check of fail rate for only flag_exp junctions! */

data junc_exp_flag_only;
   set jnc_exp_only;
   if flag_all_on=1 then output;
   keep fusion_id;
run;

proc sort data=junc_exp_flag_only nodup;
   by fusion_id; 
run;

proc sort data=flag_resids_logapn;
   by fusion_id;
run;

data exp_junc_w_resid_flag;
   merge junc_exp_flag_only (in=in1) flag_resids_logapn (in=in2);
   by fusion_id;
   if in1 and in2 then output;
run;

proc freq data=exp_junc_w_resid_flag;
   tables flag_fail_norm;
run;

/* with only flag_all_on junctions:
0 = 37504 (93.78%)
1 = 1824 (6.22%) */

/* make permenant */

data mysas.exp_junc_w_resid_flag_&sysparm.;
set exp_junc_w_resid_flag;
run;


