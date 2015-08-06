/********************************************************************************
* Need to check the distribution of the data after normalization.
* 
* libname CEGS '!MCLAB/cegs_sergey/sas_data';
* libname CEGLOCAL '!SASLOC1/cegs_sergey/sasdata';
********************************************************************************/

/* Mated */
data mated;
    retain sample_id fusion_id log_uq_apn;
    set CEGLOCAL.norm_basic_stats_m;
    sample_id = trim(line) || '_' || trim(mating_status) || trim(rep);
    keep sample_id fusion_id log_uq_apn;
    run;

proc export data=mated outfile='/tmp/mated.csv' dbms=csv replace;
putnames=yes;
run;

data _null_;
    call system('Rscript $MCLAB/cegs_sergey/r_programs/plot_normalization_distributions_v2.R /tmp/mated.csv Mated');
    run;

/* Virgin */
data virgin;
    retain sample_id fusion_id log_uq_apn;
    set CEGLOCAL.norm_basic_stats_v;
    sample_id = trim(line) || '_' || trim(mating_status) || trim(rep);
    keep sample_id fusion_id log_uq_apn;
    run;

proc export data=mated outfile='/tmp/virgin.csv' dbms=csv replace;
putnames=yes;
run;

data _null_;
    call system('Rscript $MCLAB/cegs_sergey/r_programs/plot_normalization_distributions_v2.R /tmp/virgin.csv Virgin');
    run;

/* Clean up */
proc datasets nolit;
    delete mated;
    delete virgin;
    run; quit;
