/********************************************************************************
* Want to create a clean dataset that I can send to USC.
* 
* libname cegs '!MCLAB/cegs_sergey/sas_data';
* libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
* filename mymacros '!MCLAB/cegs_sergey/sas_programs/macros';
* options SASAUTOS=(sasautos mymacros);
********************************************************************************/

/* Create design file for iterating over */
    data design_file;
        set CEGLOCAL.ccfus_norm_centered;
        keep line mating_status rep;
        run;

    proc sort data=design_file nodupkey;
        by line mating_status rep;
        run;

/* Pull out only the columns I will need so things go quicker in the macro */
    data alldata;
        set CEGLOCAL.ccfus_norm_centered;
        keep fusion_id line mating_status rep log_uq_apn uq_log_uq_center;
        run;

/* Iterate through all of the samples and prep for merging side-by-side */
    %macro make_sbs(line, mating_status, rep);

        data center_&line._&mating_status.&rep.;
            set alldata;
            if line eq "&line" and mating_status eq "&mating_status" and rep eq "&rep";
            rename uq_log_uq_center = &line._&mating_status.&rep._center;
            drop line mating_status rep log_uq_apn;
            run;

        proc sort data = center_&line._&mating_status.&rep.;
            by fusion_id;
            run;

        data norm_&line._&mating_status.&rep.;
            set alldata;
            if line eq "&line" and mating_status eq "&mating_status" and rep eq "&rep";
            rename log_uq_apn = &line._&mating_status.&rep._norm;
            drop line mating_status rep uq_log_uq_center;
            run;

        proc sort data = norm_&line._&mating_status.&rep.;
            by fusion_id;
            run;
    %mend;

    %iterdataset(dataset= design_file, function=%nrstr(%make_sbs(&line, &mating_status, &rep);));

/* Merge datasets into a side-by-side format */
    options missing = '0';
    data CEGLOCAL.norm_sbs;
        merge norm_:;
        by fusion_id;
        run;

    data CEGLOCAL.center_sbs;
        merge center_:;
        by fusion_id;
        run;

    options missing = ' ';
    data CEGLOCAL.norm_sbs;
        set merged_norm;
        run;

/* Clean up */
    proc datasets nolist;
        delete center_:;
        delete norm_:;
        delete design_file;
        delete alldata;
        run; quit;


