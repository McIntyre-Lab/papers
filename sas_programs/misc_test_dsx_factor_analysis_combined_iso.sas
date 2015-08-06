/********************************************************************************
* Perform a simple factor analysis using the collapsed isoforms.
* Looking for factor loadings as well as making sure all of the
* matrices perform well with the inversions and rotations.
********************************************************************************/

/* Combine paternal and maternal RIL ids into a single sample id */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_dsx_sbs_combine_sym;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* Perform Factor Analisys */
    proc sort data=w_sample nodupkey;
        by sample;
        run;

    options ls=250 ps=60;
    proc factor data=w_sample simple method=PRINCIPAL 
        scree rotate=varimax nfact=25 round flag=.4 outstat=factor;
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
        factor_max = max(Factor1, Factor2, Factor3, Factor4, Factor5, Factor6, Factor7, Factor8, Factor9, Factor10, Factor11, Factor12, Factor13, Factor14, Factor15, Factor16, Factor17, Factor18, Factor19, Factor20, Factor21, Factor22, Factor23, Factor24, Factor25, Factor26, Factor27, Factor28, Factor29, Factor30, Factor31);
        if factor_max lt 0.4 then flag_factor = 0;
        else if Factor1 eq factor_max then flag_factor = 1;
        else if Factor2 eq factor_max then flag_factor = 2;
        else if Factor3 eq factor_max then flag_factor = 3;
        else if Factor4 eq factor_max then flag_factor = 4;
        else if Factor5 eq factor_max then flag_factor = 5;
        else if Factor6 eq factor_max then flag_factor = 6;
        else if Factor7 eq factor_max then flag_factor = 7;
        else if Factor8 eq factor_max then flag_factor = 8;
        else if Factor9 eq factor_max then flag_factor = 9;
        else if Factor10 eq factor_max then flag_factor = 10;
        else if Factor11 eq factor_max then flag_factor = 11;
        else if Factor12 eq factor_max then flag_factor = 12;
        else if Factor13 eq factor_max then flag_factor = 13;
        else if Factor14 eq factor_max then flag_factor = 14;
        else if Factor15 eq factor_max then flag_factor = 15;
        else if Factor16 eq factor_max then flag_factor = 16;
        else if Factor17 eq factor_max then flag_factor = 17;
        else if Factor18 eq factor_max then flag_factor = 18;
        else if Factor19 eq factor_max then flag_factor = 19;
        else if Factor20 eq factor_max then flag_factor = 20;
        else if Factor21 eq factor_max then flag_factor = 21;
        else if Factor22 eq factor_max then flag_factor = 22;
        else if Factor23 eq factor_max then flag_factor = 23;
        else if Factor24 eq factor_max then flag_factor = 24;
        else if Factor25 eq factor_max then flag_factor = 25;
        else if Factor26 eq factor_max then flag_factor = 26;
        else if Factor27 eq factor_max then flag_factor = 27;
        else if Factor28 eq factor_max then flag_factor = 28;
        else if Factor29 eq factor_max then flag_factor = 29;
        else if Factor30 eq factor_max then flag_factor = 30;
        else if Factor31 eq factor_max then flag_factor = 31;
        drop factor_max Factor:;
        run;

   proc export data=flag_factor outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_dsx_factor_analysis_combined_iso.csv' dbms=csv replace;
       putnames=yes;
       run;

/* Clean up */
proc datasets nolist;
    delete factor;
    delete factor2;
    delete factor_flip;
    delete flag_factor;
    delete w_sample;
    run; quit;

