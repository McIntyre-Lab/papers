/********************************************************************************
* I am going to create a lot of different plots throughout the normalization
* process. I need to get a better idea what is going on. We are particularly
* concerned about line effects, so they will be the focus of the plots.
********************************************************************************/

/* Mated */
data mated;
    retain sample_id fusion_id uq_log_uq_center;
    set CEGLOCAL.ccfus_norm_centered;
    where mating_status eq 'M';
    sample_id = trim(line) || '_' || trim(mating_status) || trim(rep);
    keep sample_id fusion_id uq_log_uq_center;
    run;

proc export data=mated outfile='/tmp/mated.csv' dbms=csv replace;
    putnames=yes;
    run;

data _null_;
    call system(' Rscript $MCLAB/cegs_sergey/r_programs/plot_centered_distributions_v2.R /tmp/mated.csv Mated');
    run;

/* Virgin */
data virgin;
    retain sample_id fusion_id uq_log_uq_center;
    set CEGLOCAL.ccfus_norm_centered;
    where mating_status eq 'V';
    sample_id = trim(line) || '_' || trim(mating_status) || trim(rep);
    keep sample_id fusion_id uq_log_uq_center;
    run;

proc export data=virgin outfile='/tmp/virgin.csv' dbms=csv replace;
    putnames=yes;
    run;

data _null_;
    call system(' Rscript $MCLAB/cegs_sergey/r_programs/plot_centered_distributions_v2.R /tmp/virgin.csv Virgin');
    run;

/* clean up */
proc datasets nolist;
    delete mated;
    delete virgin;
    run;
