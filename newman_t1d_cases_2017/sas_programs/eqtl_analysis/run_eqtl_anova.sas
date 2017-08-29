/* Running ANOVA on IR events */

/* set libraries */

libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data';

/* import data */

proc import datafile="/scratch/lfs/patcon/jnewman/eqtls/original_data/gene_&sysparm..csv"
    out=eqtl_data_&sysparm.
    dbms=csv
    replace;
    guessingrows=827836;
run;

data eqtl_data_2_&sysparm.;
    set eqtl_data_&sysparm.;
    *if snp_id='rs72819174' then delete;
    *if event_id='SE2_0393147_IR' then delete;
    if snp_id='rs72819174' and event_id='SE2_0393147_IR' then delete;
run;

proc sort data=eqtl_data_2_&sysparm.;
    by event_id snp_id cell_type;
run;


ods listing close;
proc glimmix data=eqtl_data_2_&sysparm.;
   by event_id snp_id cell_type;
   where genotype ne .;
   class genotype sex pool subject_id;
   model log_depth = genotype | sex pool V3 / ddfm=kr;
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
  by event_id snp_id cell_type;
  var Resid;
  output out = normtest_&sysparm. probn=pnorm;
  run;

data flag_resids_&sysparm.;
  set normtest_&sysparm.;
  if pnorm = . then flag_fail_norm = .;
        else if pnorm le 0.05 then flag_fail_norm = 1;
        else flag_fail_norm = 0;
  run;

/* format P-values so SAS will import them correctly */

data anova_2_&sysparm.;
   set anova_&sysparm.;
   format ProbF best32.;
run;


/* Output CSV */

proc export data=anova_2_&sysparm.
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/anova_gene_&sysparm..csv"
    dbms=csv
    replace;
run;


proc export data=lsmean_&sysparm.
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/lsmeans_gene_&sysparm..csv"
    dbms=csv
    replace;
run;


proc export data=diffs_&sysparm.
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/diffs_gene_&sysparm..csv"
    dbms=csv
    replace;
run;

proc export data=flag_resids_&sysparm.
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/resids_gene_&sysparm..csv"
    dbms=csv
    replace;
run;
