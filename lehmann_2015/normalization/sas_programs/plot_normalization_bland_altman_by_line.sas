data tmp;
    set CEGLOCAL.ccfus_norm_centered;
    log_rpkm = log2(rpkm + 1);
    log_apn = log2(apn + 1);
    drop sample uq_ff flag_raleigh mean_log_uq_apn median_log_uq_apn uq_log_uq_apn mean_log_uq_center median_log_uq_center;
    run;

proc sort data=tmp;
    by fusion_id line mating_status rep;
    run;

%macro loop_ms(MV);

    data tmp2;
        set tmp;
        if mating_status eq &MV;
        run;

    proc means data=tmp2 noprint;
        by fusion_id line;
        output out=means mean(apn)=apn mean(rpkm)=rpkm mean(uq_apn)=uq_apn mean(log_uq_apn)=log_uq_apn mean(uq_log_uq_center)=uq_log_uq_center mean(log_rpkm)=log_rpkm mean(log_apn)=log_apn;
        run;

    proc sort data=means;
        by fusion_id line;
        run;

    proc transpose data=means out=flip;
        by fusion_id;
        var apn rpkm uq_apn log_uq_apn uq_log_uq_center log_rpkm log_apn;
        id line;
        run;

    %macro loop_type(type);

        data tmp3;
            set flip;
            where _name_ eq &type;
            rename _name_ = name;
            label _name_ = ' ';
            mating_status = &MV;
            run;

        proc export data=tmp3 outfile='/home/jfear/tmp.csv' dbms=csv replace;
        putnames = yes;
        run;

        data _null_;
            call system ('Rscript $MCLAB/cegs_sergey/r_programs/plot_normalization_bland_altman_by_line.R');
            run;

    %mend;
    %loop_type('log_apn');
    %loop_type('log_rpkm');
    %loop_type('log_uq_apn');
    %loop_type('uq_log_uq_center');
%mend;
%loop_ms('M');
%loop_ms('V');


/* Clean Up */
proc datasets nolist;
    delete tmp;
    delete tmp2;
    delete tmp3;
    delete flip;
    delete lines;
    run; quit;

