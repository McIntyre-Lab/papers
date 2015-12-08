/* import libraries */
libname con '/home/jrbnewman/concannon/sas_data';

proc sort data=con.fusions_on_gt_apn0;
   by fusion_id;
run;

proc sort data=con.fusion_q3_norm_data_all;
   by fusion_id;
run;

data fusion_q3_norm_data_all noflags_oops no_exp;
   merge con.fusion_q3_norm_data_all (in=in1) con.fusions_on_gt_apn0 (in=in2);
   by fusion_id;
   if in1 and in2 then output fusion_q3_norm_data_all;
   else if in1 then output noflags_oops;
   else output no_exp;
run;

data data_for_models;
set fusion_q3_norm_data_all;
  by fusion_id;
 
/* mhalanobis distance and low coverage*/
if name =  '2009-PC-0221' then delete; *sample 75 cd8;
if name =  '2009-PC-0144' then delete; *sample 48 cd4;
if name =  '2009-PC-0236' then delete; *sample 80;
if name =  '2009-PC-0237' then delete; *sample 80;
if name =  '2009-PC-0235' then delete; *sample 80;
 /*mahalnobis distance samples*/
  if name = '2009-PC-0101' then delete; *sample 34cd8;
if name = '2009-PC-0104' then delete; *sample35 cd8;
if name =  '2009-PC-0153' then delete; *sample 51 cd4;
if name =  '2009-PC-0200' then delete; *sample 67 cd8;
if name =  '2009-PC-0212' then delete; *sample 72 cd8;
if name = '2009-PC-0215' then delete; *sample 73  cd8;
/*funky heatmap samples*/

if name =  '2009-PC-0083' then delete; 
if name =   '2009-PC-0114' then delete;
if name =   '2009-PC-0224' then delete;
if name =   '2009-PC-0228' then delete;
run;


proc sort data=data_for_models;
   by fusion_id Name;
run;

/* Add in covariates! */

data covar_for_models;
   set con.all_covariates;
   keep name pool sex v3 cell_type subject_id;
run;

proc sort data=covar_for_models;
   by name;
run;

proc sort data=data_for_models;
by name;
run;

data data_for_models_w_cov no_exp oops;
   merge covar_for_models (in=in1) data_for_models (in=in2);
   by name;
   if in1 and in2 then output data_for_models_w_cov;
   else if in1 then output no_exp;
   else output oops;
run;


/*check testset is really all_on*/
data all_on;
  set data_for_models_w_cov;
  if flag_fusion_all_on0 = 1;
  run;

* Find order for contrasts;
ods listing;
  ods html close;

proc sort data=all_on;
  by fusion_id cell_type Name;
run;

proc freq data=all_on;
  tables cell_type;
  run;

ods listing close;
proc glimmix data=all_on;
   by fusion_id;
   class cell_type sex pool subject_id;
   model log_q3_q3_apn = cell_type | sex pool V3 / ddfm=kr;
   random pool;
   random resid / subject=subject_id;

output out=con.resid_fusions resid=resid pred=pred student=stu;


lsmeans cell_type/diff;			/*CD19		CD4		CD8*/
contrast 'CD4-CD8' cell_type		0		 1		-1	;
contrast 'CD4-CD19' cell_type		1		-1	 	 0	;
contrast 'CD8-CD19' cell_type	 	1		 0		-1	;

estimate 'CD4-CD8' cell_type		0		 1		-1	;
estimate 'CD4-CD19' cell_type		1		-1		 0	;
estimate 'CD8-CD19' cell_type		1		 0		-1	;

ods output tests3=con.anova_fusions
	lsmeans = con.lsmean_fusions
	diffs = con.diffs_fusions
	contrasts = con.contrasts_fusions
	estimates = con.estimates_fusions
;

run;
quit;

* Flag residuals;
proc univariate data = con.resid_fusions normal noprint;
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

proc freq data = flag_resids noprint;
  tables flag_fail_norm / out=resid_fail_freq;
  run;

/* Make permenant */

data con.fusions_resid_normtest;
   set flag_resids;
run;

