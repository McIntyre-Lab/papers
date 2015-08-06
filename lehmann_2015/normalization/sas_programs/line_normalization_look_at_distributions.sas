/********************************************************************************
* Need to check the distribution of the data after normalization.
* 
* libname CEGS '!MCLAB/cegs_sergey/sas_data';
* libname CEGLOCAL '!SASLOC1/cegs_sergey/sasdata';
********************************************************************************/
libname NORM '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';

/* Mated */
data mated;
    retain sample_id fusion_id log_uq_apn;
    set NORM.line_norm_basic_stats_m;
    sample_id = trim(line) || '_' || trim(mating_status);
    keep sample_id fusion_id log_uq_apn;
    run;

proc export data=mated outfile='/tmp/mated.csv' dbms=csv replace;
putnames=yes;
run;

/* Virgin */
data virgin;
    retain sample_id fusion_id log_uq_apn;
    set NORM.line_norm_basic_stats_v;
    sample_id = trim(line) || '_' || trim(mating_status);
    keep sample_id fusion_id log_uq_apn;
    run;

proc export data=mated outfile='/tmp/virgin.csv' dbms=csv replace;
putnames=yes;
run;

data _null_;
    call system('Rscript $MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/r_programs/plot_normalization_distributions_line.R /tmp/mated.csv /tmp/virgin.csv');
    run;

/* Clean up */
proc datasets nolist;
    delete mated;
    delete virgin;
    run; quit;
