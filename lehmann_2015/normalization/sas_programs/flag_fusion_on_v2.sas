/* Remove Samples that have low expression */
data alldata;
    set CEGLOCAL.ccfus_stack;
    run;

proc sort data=alldata;
    by line mating_status rep;
    run;

proc sort data=CEGS.flag_sample_on;
    by line mating_status rep;
    run;

proc sort data=CEGS.flag_contaminated;
    by line mating_status rep;
    run;

data clean;
    merge alldata (in=in1) CEGS.flag_sample_on (in=in2) CEGS.flag_contaminated (in=in3);
    by line mating_status rep;
    if flag_sample_on = 1 and flag_contaminated = 0;
    run;
    
/* Flag Fusions with any expression */
data expressed;
    set clean;
    if APN > 0 then flag_expressed = 1; else flag_expressed =0;
    keep fusion_id line mating_status rep flag_expressed;
    run;

/* Create Dataset with counts of detected exons */
    proc freq data=expressed noprint;
        by line mating_status rep;
        table flag_expressed / out=freqs;
        run;

    data CEGS.detected_exon_cnts;
        set freqs;
        where flag_expressed = 1;
        rename count = detected_exon_cnts;
        label count = ' ';
        drop flag_expressed percent;
        run;

/* Summarize to fusion*line*mv level and flag a fusion as on if it is in >50% of line*reps */
proc sort data=expressed;
    by fusion_id line mating_status;
    run;

proc means data=expressed noprint;
    by fusion_id line mating_status;
    output out=means mean(flag_expressed)=mean_flag_expressed;
    run;

data CEGS.flag_fusion_on;
    set means;
    if mean_flag_expressed >= .5 then flag_fusion_on = 1; else flag_fusion_on = 0;
    keep fusion_id line mating_status flag_fusion_on;
    run;

/* Flag a Fusion to be dropped if not on in 90% of of line */
proc sort data=CEGS.flag_fusion_on;
    by fusion_id mating_status;
    run;

proc means data=CEGS.flag_fusion_on noprint;
    by fusion_id mating_status;
    output out=sums sum(flag_fusion_on)=sum_flag_fusion_on;
    run;

data CEGS.flag_drop_fusion;
    set sums;
    prop = sum_flag_fusion_on / _freq_;
    if prop < .90 then flag_drop_fusion = 1; else flag_drop_fusion = 0;
    keep fusion_id mating_status flag_drop_fusion;
    run;

proc sort data=CEGS.flag_drop_fusion;
    by mating_status;
    run;

proc freq data=CEGS.flag_drop_fusion;
    by mating_status;
    table flag_drop_fusion;
    run;
    
proc datasets nolist;
    delete alldata;
    delete clean;
    delete means;
    delete sums;
    delete expressed;
    run; quit;
