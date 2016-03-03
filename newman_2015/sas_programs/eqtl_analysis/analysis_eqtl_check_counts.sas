/* Output frequency of each combination of gene, feature, and snp */

libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';

proc sort data=eqtl.eqtl_data_for_models_&sysparm.;
     by gene_id feature_id snp_id cell_type;
run;

proc freq data=eqtl.eqtl_data_for_models_&sysparm. noprint;
    by gene_id feature_id snp_id;
    tables cell_type / out=eqtl_data_&sysparm.;
    run;
    
    
proc export data=eqtl_data_&sysparm.
    outfile="/scratch/lfs/patcon/jnewman/sas_analysis/eqtl_data_check_&sysparm..csv"
    dbms=csv replace;
    run;

