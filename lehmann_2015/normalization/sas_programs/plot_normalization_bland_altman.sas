data tmp;
    set CEGLOCAL.ccfus_norm_centered;
    log_rpkm = log2(rpkm + 1);
    log_apn = log2(apn + 1);
    drop sample uq_ff flag_raleigh mean_log_uq_apn median_log_uq_apn uq_log_uq_apn mean_log_uq_center median_log_uq_center;
    run;

proc sort data=tmp;
    by fusion_id line mating_status rep;
    run;

%macro loop_line(line);
    %macro loop_ms(MV);

        data tmp2;
            set tmp;
            if line eq &line and mating_status eq &MV;
            run;

        proc transpose data=tmp2 out=flip prefix=rep;
            by fusion_id line mating_status ;
            var apn rpkm uq_apn log_uq_apn uq_log_uq_center log_rpkm log_apn;
            id rep;
            run;

        %macro loop_type(type);

            data tmp3;
                set flip;
                where _name_ eq &type;
                rename _name_ = name;
                label _name_ = ' ';
                run;

            proc export data=tmp3 outfile='/home/jfear/tmp.csv' dbms=csv replace;
            putnames = yes;
            run;

            data _null_;
                call system ('Rscript $MCLAB/cegs_sergey/r_programs/plot_normalization_bland_altman.R');
                run;

        %mend;
        %loop_type('apn');
        %loop_type('rpkm');
        %loop_type('uq_apn');
        %loop_type('log_apn');
        %loop_type('log_rpkm');
        %loop_type('log_uq_apn');
        %loop_type('uq_log_uq_center');
    %mend;
    %loop_ms('M');
    %loop_ms('V');
%mend;

data lines;
    set CEGS.combined_design_by_rep;
    drop mating_status rep;
    run;

proc sort data=lines nodupkey;
    by line;
    run;

%iterdataset(dataset=lines,function=%nrstr(%loop_line("&line");));

/* Clean Up */
proc datasets nolist;
    delete tmp;
    delete tmp2;
    delete tmp3;
    delete flip;
    delete lines;
    run; quit;

