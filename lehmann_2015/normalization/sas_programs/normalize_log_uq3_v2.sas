/********************************************************************************
* Normalization using the log2 of the upper quartile adjusted APN.
* Normalization is performed separately for Mated and Virgin. Normalizations
* include the coverage counts on Junctions.
********************************************************************************/
/* Create clean dataset */
    * Remove samples that have little expression or are contaminated;
    proc sort data=CEGLOCAL.ccfus_stack;
        by line mating_status rep;
        run;

    proc sort data=CEGS.flag_sample_on;
        by line mating_status rep;
        run;

    proc sort data=CEGS.flag_contaminated;
        by line mating_status rep;
        run;

    data clean1;
        merge CEGLOCAL.ccfus_stack (in=in1) CEGS.flag_sample_on (in=in2) CEGS.flag_contaminated (in=in3);
        by line mating_status rep;
        if flag_sample_on = 1 and flag_contaminated = 0;
        drop flag_sample_on flag_contaminated;
        run;

    * Remove fusions with little expression;
        proc sort data=clean1;
            by fusion_id mating_status;
            run;

        proc sort data=CEGS.flag_drop_fusion;
            by fusion_id mating_status;
            run;

        data clean2;
            merge clean1 (in=in1) CEGS.flag_drop_fusion;
            by fusion_id mating_status;
            if flag_drop_fusion = 0;
            drop flag_drop_fusion;
            run;

/* Create Clean Junction dataset */
    proc sort data=CEGLOCAL.junc_cnts;
        by line mating_status rep;
        run;

    data junc_clean1;
        merge CEGLOCAL.junc_cnts (in=in1) CEGS.flag_sample_on (in=in2) CEGS.flag_contaminated (in=in3);
        by line mating_status rep;
        if flag_sample_on = 1 and flag_contaminated = 0;
        rename junc_reads_in_region = reads_in_region;
        drop flag_sample_on flag_contaminated;
        run;
        
/* Combine Fusion coverage counts and junction coverage counts */
    data combined_clean;
        set clean2 (in=in1) junc_clean1 (in=in2);
        if in2 then flag_junc = 1; else flag_junc = 0;
        run;

    proc sort data=combined_clean;
        by line mating_status rep;
        run;

/* Calculate With-in sample basic statistics */
    proc means data=combined_clean noprint;
        class line mating_status rep;
        var reads_in_region;
        output out=mapped_reads sum=sum_mapped q3=q3 median=median;
        run;

    data totals;
        set mapped_reads;
        where _type_ = 7;
        keep line mating_status rep sum_mapped q3 median;
        run;

    /* Summarize statistics to experiment level */
        proc sort data=totals;
            by mating_status ;
            run;

        proc means data=totals median ;
            by mating_status;
            var sum_mapped q3 median;
            run;

    /* Merge statistics onto dataset */
    proc sort data=clean2;
        by line mating_status rep;
        run;

    proc sort data=totals;
        by line mating_status rep;
        run;

    data stack oops;
        merge clean2 (in=in1) totals (in=in2);
        by line mating_status rep;
        if in1 and in2 then output stack;
        else output oops; * 0 obs yay!;
        run;

    /* Merge on junction junc_mapped_reads to get total mapped_reads */
        proc sort data=stack;
            by line mating_status rep;
            run;

        data junc_mapped;
            set junc_clean1;
            keep line mating_status rep junc_mapped_reads;
            run;

        proc sort data=junc_mapped nodupkey;
            by line mating_status rep;
            run;

        data stack2 oops;
            merge stack (in=in1) junc_mapped (in=in2);
            by line mating_status rep;
            total_mapped_reads = mapped_reads + junc_mapped_reads;
            if in1 then output stack2;
            else output oops; * 0 obs;
            run;

    /* NOTE MAKE SURE TO CHANGE THE UQ MULTIPLIER BY THE NEW Q3 */
    data stack3_m;
        retain fusion_id;
        set stack2 (where=(mating_status = 'M'));
        uq_apn = (apn/q3)*33.8119518;
        uq_ff = 33.8119518/q3;
        log_uq_apn = log2(UQ_apn + 1);
        run;

    data stack3_v;
        retain fusion_id;
        set stack2 (where=(mating_status = 'V'));
        uq_apn = (apn/q3)*33.8645833;
        uq_ff = 33.8645833/q3;
        log_uq_apn = log2(UQ_apn + 1);
        run;

    /* Make sure that all fusions are present */
    proc freq data=stack3_m noprint;
        table fusion_id/ out=freqs;
        run; * 31079 obs OK!!! same as freqs from previous script;

    data CEGS.norm_basic_stats_m;
        set stack3_m;
        run;

    data CEGLOCAL.norm_basic_stats_m;
        set stack3_m;
        run;

    proc freq data=stack3_v noprint;
        table fusion_id/ out=freqs;
        run; * 32722 obs OK!!! same as freqs from prevsious script;

    data CEGS.norm_basic_stats_v;
        set stack3_v;
        run;

    data CEGLOCAL.norm_basic_stats_v;
        set stack3_v;
        run;

/* Clean up */
proc datasets nolist;
    delete clean1;
    delete clean2;
    delete combined_clean;
    delete freqs;
    delete junc_clean1;
    delete junc_mapped;
    delete mapped_reads;
    delete oops;
    delete stack;
    delete stack2;
    delete stack3_m;
    delete stack3_v;
    delete totals;
    run; quit;
