/********************************************************************************
* Using the Mahalnobis distance, flag outliters and look at some box plots to
* determine if distributions look very different.
*
* libname cegs '!MCLAB/cegs_sergey/sas_data';
* libname ceglocal '!SASLOC1/cegs_sergey/sasdata';
********************************************************************************/

data mal;
    set CEGS.mahalanobis_distance;
    run;

/* Use design file to merge on rep */
    data design;
        set CEGS.combined_design_by_rep;
        _name_ = trim(line) || '_' || trim(mating_status) || trim(rep) || '_center';
        run;

    proc sort data=design;
        by _name_;
        run;

    proc sort data=mal;
        by _name_;
        run;

    data merged;
        merge design (in=in1) mal (in=in2);
        by _name_;
        if in1 and in2;
        run;

/* Create Flag_mahalanobis_outlier */
    data CEGS.flag_mahalanobis_outlier;
        set merged;
        if mahal_distances ge 4 then flag_mahalanobis_outlier = 1;
        else flag_mahalanobis_outlier = 0;
        drop _name_ mahal_distances;
        run;

/* Merge on counts and export files for plots */
    proc sort data=CEGS.flag_mahalanobis_outlier;
        by line mating_status rep;
        run;

    proc sort data= CEGLOCAL.ccfus_norm_centered;
        by line mating_status rep;
        run;

    data merged2;
        merge CEGS.flag_mahalanobis_outlier CEGLOCAL.ccfus_norm_centered;
        by line mating_status rep;
        keep line mating_status rep fusion_id uq_log_uq_center flag_mahalanobis_outlier;
        run;

    proc export data=merged2 outfile='/tmp/outlier.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Run Rscript to make plots */
    data _null_;
        call system('Rscript $MCLAB/cegs_sergey/r_programs/plot_mahalanobis_outliers.R /tmp/outlier.csv');
        run;

/* Summarize to line level */
proc freq data=CEGS.flag_mahalanobis_outlier noprint;
    by line mating_status;
    table flag_mahalanobis_outlier /out=freqs;
    run;
    
data freqs2;
    set freqs;
    where flag_mahalanobis_outlier = 1;
    run;

proc sort data=freqs2;
    by DESCENDING percent;
    run;

