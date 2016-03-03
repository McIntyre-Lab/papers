/* Stack eQTL datasets and prep for ANOVA */

libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';

/* Merge in covariate information */

data covariates;
   set eqtl.all_covariates;
   keep Name cell_type subject_id pool sex V3;
run;

proc sort data=covariates;
   by Name;
run;

proc sort data=eqtl.all_exp_data_w_snps_&sysparm.;
   by Name;
run;


/* Merge covariates with eQTL data */

data eqtl_data_w_covar;
   merge covariates (in=in1) eqtl.all_exp_data_w_snps_&sysparm. (in=in2);
   by Name;
   if in1 and in2;
run;


/* Drop low coverage samples */


data eqtl.eqtl_data_for_models_&sysparm.;
   set eqtl_data_w_covar;
   if name =  '2009-PC-0221' then delete; *sample 75 cd8;
   if name =  '2009-PC-0144' then delete; *sample 48 cd4;
   if name =  '2009-PC-0236' then delete; *sample 80;
   if name =  '2009-PC-0237' then delete; *sample 80;
   if name =  '2009-PC-0235' then delete; *sample 80;
   run;

proc sort data=eqtl.eqtl_data_for_models_&sysparm. nodup;
    by feature_id snp_id cell_type;
run;

