/********************************************************************************
* Perform a simple factor analysis using the collapsed isoforms.
* Looking for factor loadings as well as making sure all of the
* matrices perform well with the inversions and rotations.
********************************************************************************/

/* Combine paternal and maternal RIL ids into a single sample id */
    data w_sample;
        retain sample;
        format sample $ 20.;
        set SEM.dsrp_sex_det_sbs_gene_level_sym;
        sample = strip(patRIL) || '_' || strip(matRIL);
        drop patRIL matRIL;
        run;

/* Perform Factor Analisys */
    proc sort data=w_sample nodupkey;
        by sample;
        run;

    options ls=250 ps=60;
    proc factor data=w_sample simple method=PRINCIPAL 
        scree rotate=varimax nfact=8 round flag=.4 outstat=factor;
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
        factor_max = max(Factor1,Factor2,Factor3,Factor4,Factor5, Factor6, Factor7, Factor8);
        if factor_max lt 0.4 then flag_factor = 0;
        else if Factor1 eq factor_max then flag_factor = 1;
        else if Factor2 eq factor_max then flag_factor = 2;
        else if Factor3 eq factor_max then flag_factor = 3;
        else if Factor4 eq factor_max then flag_factor = 4;
        else if Factor5 eq factor_max then flag_factor = 5;
        else if Factor6  eq factor_max then flag_factor = 6 ;
        else if Factor7  eq factor_max then flag_factor = 7 ;
        else if Factor8  eq factor_max then flag_factor = 8 ;
        drop factor_max Factor:;
        run;

    proc export data=flag_factor outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_genes.csv' dbms=csv replace;
        putnames=yes;
        run;


/* Create Network Plots In Python */
data _null_;
    call system('python $MCLAB/cegs_sem_sd_paper/scripts/plot_fa_grid.py --factor $MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_genes.csv -o $MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_genes.png -fs 14 -cs 20');
    run;

/* Clean up */
proc datasets;
    delete factor;
    delete factor2;
    delete factor_flip;
    delete flag_factor;
    delete w_sample;
    run; quit;

