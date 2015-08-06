/********************************************************************************
* I am going to create a lot of different plots throughout the normalization
* process. I need to get a better idea what is going on. We are particularly
* concerned about line effects, so they will be the focus of the plots.
********************************************************************************/
/* Pull in normalized data for Mated and Virgin */
    data alldata;
        set CEGLOCAL.ccfus_norm_stack_m CEGLOCAL.ccfus_norm_stack_v;
        keep fusion_id line mating_status sample rep apn rpkm uq_apn log_uq_apn uq_ff;
        run;

/* Merge on flag_raleigh */
    proc sort data=alldata;
        by line;
        run;

    proc sort data=CEGS.flag_raleigh;
        by line;
        run;

    data alldata2;
        merge alldata (in=in1) CEGS.flag_raleigh (in=in2);
        by line;
        if in1;
        run;

/* Calculate measures of central tendancy */
    proc sort data=alldata;
        by line mating_status rep;
        run;

    proc means data=alldata noprint;
        by line mating_status rep;
        output out=means 
        mean(log_uq_apn)=mean_log_uq_apn 
        median(log_uq_apn)=median_log_uq_apn 
        q3(log_uq_apn)=uq_log_uq_apn 
        q3(rpkm)=uq_rpkm 
        ;
        run;

    proc sort data=means;
        by line mating_status rep;
        run;

    proc sort data=alldata2;
        by line mating_status rep;
        run;

    data alldata3;
        merge alldata2 (in=in1) means (in=in2);
        by line mating_status rep;
        if in1;
        drop _type_ _freq_;
        run;

/* Calculate Centered Values */
    data CEGLOCAL.ccfus_norm_centered;
        set alldata3;
        mean_log_uq_center = log_uq_apn - mean_log_uq_apn;
        median_log_uq_center = log_uq_apn - median_log_uq_apn;
        uq_log_uq_center = log_uq_apn - uq_log_uq_apn;
        uq_rpkm_center = rpkm - uq_rpkm;
        run;

    data CEGS.ccfus_norm_centered;
        set CEGLOCAL.ccfus_norm_centered;
        run;

/* Clean Up */
proc datasets nolist;
    delete alldata;
    delete alldata2;
    delete alldata3;
    delete means;
    run; quit;
