/* ANOVA model. APN>0

*/

libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';

libname ribo 'S:\McIntyre_Lab\arbeitman\arbeitman_ribotag\sas_data';
libname dmel 'S:\McIntyre_Lab\useful_dmel_data\flybase551\sasdata';

proc sort data=ribo.all_coverage_counts_with_key;
  by fusion_id;
  run;

proc sort data=ribo.on_calls_gt_apn0;
  by fusion_id;
  run;

data ribo_all;
  merge ribo.all_coverage_counts_with_key (in=in1) ribo.on_calls_gt_apn0 (in=in2);
  by fusion_id;
  if in1;
  run;

proc freq data = ribo_all;
  tables flag_fusion_all_on0;
  run; *flag=1 37.7%;

data ribo.ribo_all;
  set ribo_all;
  run;

data ribo_fusions_on;
  set ribo_all;
  if flag_fusion_all_on0=1;
  run;

/*Determine order for contrasts*/
proc freq data = ribo_fusions_on;
  tables trt ;
  run;

/*Contrasts of interest for mixed model*/
ods listing close;
proc mixed data = ribo_fusions_on;
  by fusion_id;
  class trt;
  model logRPKM= trt / htype=1 outp=resid;

lsmeans trt/diff;			/*IPfemale	IPmale		InputFemale	InputMale*/
contrast 'IPmale-IPfemale' trt		1		-1		0		0		;
contrast 'IPmale-InputMale' trt		0		1		0		-1		;
contrast 'IPfemale-InputFemale' trt	1		0		-1		0		;
contrast 'InputMale-InputFemale' trt	0		0		1		-1		;

estimate 'IPmale-IPfemale' trt		1		-1		0		0		;
estimate 'IPmale-InputMale' trt		0		1		0		-1		;
estimate 'IPfemale-InputFemale' trt	1		0		-1		0		;
estimate 'InputMale-InputFemale' trt	0		0		1		-1		;

ods output tests1 = anova
	lsmeans = lsmean
	diffs = diffs
	contrasts = contrasts
	estimates = estimates;
run;
quit;
ods listing;

/*Make datasets permanent*/
data ribo.output_anova_apn0 ;
set anova ;
run ;
data ribo.output_lsmean_apn0 ;
set lsmean ;
run ;
data ribo.output_diffs_apn0 ;
set diffs ;
run ;
data ribo.output_contrasts_apn0 ;
set contrasts ;
run ;
data ribo.output_estimates_apn0 ;
set estimates ;
run ;


/* Flag residuals that fail normality */

proc univariate data = resid normal noprint;
  by fusion_id;
  var Resid;
  output out = normtest probn=pnorm;
  run;

data flag_resids;
  set normtest;
  if pnorm = . then flag_fail_norm = .;
  	else if pnorm le 0.05 then flag_fail_norm = 1;
	else flag_fail_norm = 0;
  run;

proc freq data = flag_resids;
  tables flag_fail_norm;
  run; *46.85% failed normality;

/*Make dataset permanent*/
data ribo.output_flag_resids_apn0;
  set flag_resids;
  run;



