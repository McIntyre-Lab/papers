/********************************************************************************
* I am going to create a lot of different plots throughout the normalization
* process. I need to get a better idea what is going on. We are particularly
* concerned about line effects, so they will be the focus of the plots.
********************************************************************************/
libname NORM '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';

/* Mated */
data mated;
    retain sample_id fusion_id uq_log_uq_center;
    set NORM.line_norm_centered;
    where mating_status eq 'M';
    sample_id = trim(line) || '_' || trim(mating_status);
    keep sample_id fusion_id uq_log_uq_center;
    run;

proc export data=mated outfile='/tmp/mated.csv' dbms=csv replace;
    putnames=yes;
    run;


/* Virgin */
data virgin;
    retain sample_id fusion_id uq_log_uq_center;
    set NORM.line_norm_centered;
    where mating_status eq 'V';
    sample_id = trim(line) || '_' || trim(mating_status);
    keep sample_id fusion_id uq_log_uq_center;
    run;

proc export data=virgin outfile='/tmp/virgin.csv' dbms=csv replace;
    putnames=yes;
    run;

data _null_;
    call system(' Rscript $MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/r_programs/plot_centered_distributions_line.R /tmp/mated.csv /tmp/virgin.csv ');
    run;

/* clean up */
proc datasets nolist;
    delete mated;
    delete virgin;
    run;
