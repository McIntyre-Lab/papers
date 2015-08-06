/********************************************************************************
* Perform a simple factor analysis using all of the genes in the sex det
* pathway.
********************************************************************************/

/* Combine paternal and maternal RIL ids into a single sample id */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_sex_det_sbs_symbol ;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* Perform Factor Analisys */
    proc sort data=w_sample nodupkey;
    by sample;
    run;

    options ls=250 ps=60;
    proc factor data=w_sample simple method=PRINCIPAL 
    scree rotate=varimax nfact=17 round flag=.4 outstat=factor;
    var : ;
    run;

/* Pull out factors and create a summary dataset */
    data factor2;
        set factor;
        if _type_ = 'PATTERN';
        run;

    proc transpose data=factor2 out=factor_flip;
    var :;
    id _name_;
    run;

    data flag_factor;
        set factor_flip;
        rename _name_ = symbol;
        factor_max = max(Factor1,Factor2,Factor3,Factor4,Factor5, Factor6, Factor7, Factor8, Factor9, Factor10, Factor11, Factor12, Factor13, Factor14, Factor15, Factor16, Factor17);
        if factor_max lt 0.4 then flag_factor = 0;
        else if Factor1 eq factor_max then flag_factor = 1;
        else if Factor2 eq factor_max then flag_factor = 2;
        else if Factor3 eq factor_max then flag_factor = 3;
        else if Factor4 eq factor_max then flag_factor = 4;
        else if Factor5 eq factor_max then flag_factor = 5;
        else if Factor6  eq factor_max then flag_factor = 6 ;
        else if Factor7  eq factor_max then flag_factor = 7 ;
        else if Factor8  eq factor_max then flag_factor = 8 ;
        else if Factor9  eq factor_max then flag_factor = 9 ;
        else if Factor10 eq factor_max then flag_factor = 10;
        else if Factor11 eq factor_max then flag_factor = 11;
        else if Factor12 eq factor_max then flag_factor = 12;
        else if Factor13 eq factor_max then flag_factor = 13;
        else if Factor14 eq factor_max then flag_factor = 14;
        else if Factor15 eq factor_max then flag_factor = 15;
        else if Factor16 eq factor_max then flag_factor = 16;
        else if Factor17 eq factor_max then flag_factor = 17;
        drop factor_max Factor:;
        run;

    proc export data=flag_factor outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_all_iso.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Create Network Plots In Python NOT DONE YET */
    /*
        data _null_;
            call system('python $MCLAB/cegs_sem_sd_paper/scripts/graph_network_factor_analysis_uq_log_uq_apn.py');
            run;
    */
