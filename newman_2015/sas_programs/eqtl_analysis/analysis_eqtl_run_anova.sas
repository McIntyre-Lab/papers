/* Run ANOVA */

libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';

data eqtl_data_for_models_&sysparm.;
   set eqtl.eqtl_data_for_models_&sysparm.;
   drop gene_id;
run;

proc sort data=eqtl_data_for_models_&sysparm. nodup;
    by feature_id snp_id cell_type Name;
run;

ods listing close;
proc glimmix data=eqtl_data_for_models_&sysparm.;
   by feature_id snp_id cell_type;
   where genotype ne .;
   class genotype sex pool subject_id;
   model log_measurement = genotype | sex pool V3 / ddfm=kr;
   random pool;
   random resid / subject=subject_id;


output out=resid_&sysparm. resid=resid pred=pred student=stu;
lsmeans genotype/diff;

ods output tests3=anova_&sysparm.
        lsmeans = lsmean_&sysparm.
        diffs = diffs_&sysparm.;
run;
quit;

 * Flag residuals;
 proc univariate data = resid_&sysparm. normal noprint;
   by feature_id snp_id cell_type;
   var Resid;
   output out = normtest_&sysparm. probn=pnorm;
   run;
 
 data flag_resids_&sysparm.;
   set normtest_&sysparm.;
   if pnorm = . then flag_fail_norm = .;
         else if pnorm le 0.05 then flag_fail_norm = 1;
         else flag_fail_norm = 0;
   run;
 

/* Make permenant */

data eqtl.anova_&sysparm. ;
    set anova_&sysparm.;
run;

data eqtl.lsmean_&sysparm.;
    set lsmean_&sysparm.;
run;

data eqtl.diffs_&sysparm.;
    set diffs_&sysparm.;
run;


data eqtl.resid_&sysparm.;
    set normtest_&sysparm.;
run;



