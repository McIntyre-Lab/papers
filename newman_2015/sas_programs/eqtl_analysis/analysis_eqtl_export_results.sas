
/* Stack ANOVA outputs and export */

libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';


data eqtl.anova_results;
   set eqtl.anova_: ;
run;

data eqtl.lsmeans_results;
   set eqtl.lsmeans_: ;
run;

data eqtl.diffs_results;
   set eqtl.diffs_: ;
run;

data eqtl.resid_results;
   set eqtl.resid_: ;
run;



 /* format P-values so SAS will import them correctly */
 
 data anova_results;
    set eqtl.anova_results;
    format ProbF best32.;
 run;






/* Output CSV */

proc export data=anova_results
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/eqtl_anova.csv"
    dbms=csv
    replace;
run;

proc export data=eqtl.lsmeans_results
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/eqtl_lsmeans.csv"
    dbms=csv
    replace;
run;

proc export data=eqtl.diffs_results
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/eqtl_diffs.csv"
    dbms=csv
    replace;
run;

proc export data=eqtl.resid_results
    outfile="/scratch/lfs/patcon/jnewman/eqtls/results/eqtl_resids.csv"
    dbms=csv
    replace;
run;


