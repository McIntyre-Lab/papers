/********************************************************************************
* I am going to create a lot of different plots throughout the normalization
* process. I need to get a better idea what is going on. We are particularly
* concerned about line effects, so they will be the focus of the plots.
********************************************************************************/
/* Pull in normalized data for Mated and Virgin */
    data alldata;
        set CEGLOCAL.line_norm_basic_stats_m CEGLOCAL.line_norm_basic_stats_v;
        keep fusion_id line mating_status apn uq_apn log_uq_apn uq_ff;
        run;

/* Calculate measures of central tendancy */
    proc sort data=alldata;
        by line mating_status ;
        run;

    proc means data=alldata noprint;
        by line mating_status ;
        output out=means 
        mean(log_uq_apn)=mean_log_uq_apn 
        median(log_uq_apn)=median_log_uq_apn 
        q3(log_uq_apn)=uq_log_uq_apn 
        ;
        run;

    proc sort data=means;
        by line mating_status ;
        run;

    proc sort data=alldata;
        by line mating_status ;
        run;

    data alldata2;
        merge alldata (in=in1) means (in=in2);
        by line mating_status ;
        if in1;
        drop _type_ _freq_;
        run;

/* Calculate Centered Values */
    data CEGLOCAL.line_norm_centered;
        set alldata2;
        mean_log_uq_center = log_uq_apn - mean_log_uq_apn;
        median_log_uq_center = log_uq_apn - median_log_uq_apn;
        uq_log_uq_center = log_uq_apn - uq_log_uq_apn;
        run;

    data CEGS.line_norm_centered;
        set CEGLOCAL.line_norm_centered;
        run;

/* Clean Up */
proc datasets nolist;
    delete alldata;
    delete alldata2;
    delete means;
    run; quit;
