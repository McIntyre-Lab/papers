/* Remove Samples that have low expression */
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

    data clean;
        merge alldata (in=in1) CEGS.flag_sample_onk (in=in2) CEGS.flag_contaminatedk (in=in3);
        by line mating_status rep;
        if flag_sample_on = 1 and flag_contaminated = 0;
        keep line mating_status rep fusion_id read_length region_length region_depth;
        run;
    
/* Summarize to line*mating_status */
    proc sort data=clean;
        by fusion_id line mating_status;
        run;
        
    proc means data=clean noprint;
        by fusion_id line mating_status;
        output out=sums sum(region_depth)=sum_region_depth;
        run;

/* Merge on fusion annotation and calculate line*mv level apn */
    data annotation;
        set clean;
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
        if apn > 0 then flag_expressed = 1; else flag_expressed =0;
        drop read_length region_length sum_region_depth;
        run;

/* Flag a Fusion to be dropped if not on in 90% of of line */
    proc sort data=apn;
        by fusion_id mating_status;
        run;

    proc freq data=apn noprint;
        by fusion_id mating_status;
        table flag_expressed /out=freqs;
        run;

    data CEGS.flag_drop_fusion_by_line;
        set freqs;
        where flag_expressed = 1;
        if percent < 90 then flag_drop_fusion = 1; else flag_drop_fusion = 0;
        keep fusion_id mating_status flag_drop_fusion;
        run;

    proc sort data=CEGS.flag_drop_fusion_by_line;
        by mating_status;
        run;

    proc freq data=CEGS.flag_drop_fusion_by_line;
        by mating_status;
        table flag_drop_fusion;
        run;
    
/* Clean up */
proc datasets nolist;
    delete alldata;
    delete annotation;
    delete apn;
    delete clean;
    delete freqs;
    delete merged;
    delete sums;
    run; quit;
