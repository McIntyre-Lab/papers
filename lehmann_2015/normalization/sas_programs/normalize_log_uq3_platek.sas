/********************************************************************************
* Normalization using the log2 of the upper quartile adjusted APN.
* Normalization is performed separately for Mated and Virgin. Normalizations
* include the coverage counts on Junctions.
********************************************************************************/
/* Create clean dataset */
    * Remove samples that have little expression or are contaminated;
    data alldata;
        set CEGLOCAL.platek_ccfus_stack CEGLOCAL.ccfus_stack;
        run;

    proc sort data=alldata;
        by line mating_status rep;
        run;

    proc sort data=CEGS.flag_sample_onk;
        by line mating_status rep;
        run;

    proc sort data=CEGS.flag_contaminatedk;
        by line mating_status rep;
        run;

    data clean1;
        merge alldata (in=in1) CEGS.flag_sample_onk (in=in2) CEGS.flag_contaminatedk (in=in3);
        by line mating_status rep;
        if flag_sample_on = 1 and flag_contaminated = 0;
        drop flag_sample_on flag_contaminated;
        run;

    * Remove fusions with little expression;
        proc sort data=clean1;
            by fusion_id mating_status;
            run;

        proc sort data=CEGS.flag_drop_fusion_by_line;
            by fusion_id mating_status;
            run;

        data clean2;
            merge clean1 (in=in1) CEGS.flag_drop_fusion_by_line;
            by fusion_id mating_status;
            if flag_drop_fusion = 0;
            drop flag_drop_fusion;
            run;

    * Summarize to line*mating_status;
        proc sort data=clean2;
            by fusion_id line mating_status;
            run;
            
        proc means data=clean2 noprint;
            by fusion_id line mating_status;
            output out=sums sum(region_depth)=sum_region_depth;
            run;

    * Merge on fusion annotation and calculate line*mv level apn;
        data annotation;
            set clean2;
            keep fusion_id read_length region_length;
            run;

        proc sort data=annotation nodupkey;
            by fusion_id;
            run;

        proc sort data=sums;
            by fusion_id ;
            run;

        data merged;
            merge annotation (in=in1) sums (in=in2);
            by fusion_id ;
            drop _type_ _freq_;
            run;

        data apn;
            set merged;
            apn = sum_region_depth / region_length;
            reads_in_region = sum_region_depth / read_length;
            run;

/* Calculate With-in sample basic statistics */
    proc means data=apn noprint;
        class line mating_status;
        var reads_in_region;
        output out=mapped_reads sum=sum_mapped q3=q3 median=median;
        run;

    data totals;
        set mapped_reads;
        where _type_ = 3;
        keep line mating_status sum_mapped q3 median;
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
    proc sort data=apn;
        by line mating_status ;
        run;

    proc sort data=totals;
        by line mating_status ;
        run;

    data stack oops;
        merge apn (in=in1) totals (in=in2);
        by line mating_status;
        if in1 and in2 then output stack;
        else output oops; * 0 obs yay!;
        run;

    /* NOTE MAKE SURE TO CHANGE THE UQ MULTIPLIER BY THE NEW Q3 */
    data stack2_m;
        retain fusion_id;
        set stack (where=(mating_status = 'M'));
        uq_ff = 156.9631579/q3;
        uq_apn = apn*uq_ff;
        log_uq_apn = log2(UQ_apn + 1);
        run;

    data stack2_v;
        retain fusion_id;
        set stack (where=(mating_status = 'V'));
        uq_ff = 156.7578947/q3;
        uq_apn = apn*uq_ff;
        log_uq_apn = log2(UQ_apn + 1);
        run;

    /* Make sure that all fusions are present */
    proc freq data=stack2_m noprint;
        table fusion_id/ out=freqs;
        run; * 37422 obs OK!!! same as freqs from previous script;

    data CEGS.line_norm_basic_stats_m;
        set stack2_m;
        run;

    data CEGLOCAL.line_norm_basic_stats_m;
        set stack2_m;
        run;

    proc freq data=stack2_v noprint;
        table fusion_id/ out=freqs;
        run; * 35676 obs OK!!! same as freqs from prevsious script;

    data CEGS.line_norm_basic_stats_v;
        set stack2_v;
        run;

    data CEGLOCAL.line_norm_basic_stats_v;
        set stack2_v;
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
