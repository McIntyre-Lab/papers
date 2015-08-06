/*******************************************************************************
* Filename: cegsV_identify_best_model_adding_genes_nofilter_y2_BIC12.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: After doing some simulations, there is a TIER of 10% until we
* get with a BIC greater in magnitude than 12.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Get the baseline value */
    data _null_;
        set SEM.cegsV_ag_yp2_stack_bic ;
        where model eq 'Model_Baseline';
        if _n_ eq 1 then do;
        call symput('Baseline', BIC);
        end;
        run;

/* How many models*gene imporved fit */
    data models;
        set SEM.cegsV_ag_yp2_stack_bic ;
        if BIC le &BASELINE - 12;
        run; * 12,565 models better than baseline;

/* For each gene figure out the lowest BIC */
    proc sort data=SEM.cegsV_ag_yp2_stack_bic ;
        by gene;
        run;

    proc means data=SEM.cegsV_ag_yp2_stack_bic noprint;
        by gene;
        var BIC;
        output out=mins min=min;
        run;

    data mins2;
        set mins;
        rename min=BIC;
        keep gene min;
        run;

    * Merge lowest BIC back to model to figure out which model has lowest BIC;
    proc sort data=mins2;
        by gene BIC;
        run;

    proc sort data=SEM.cegsV_ag_yp2_stack_bic ;
        by gene BIC;
        run;

    data merged;
        merge SEM.cegsV_ag_yp2_stack_bic  (in=in1) mins2 (in=in2);
        by gene BIC;
        if in1 and in2;
        diff = &Baseline - BIC;
        if diff < 12 then model = 'Model_Baseline';
        drop diff;
        run; 

    proc sort data=merged;
        by model;
        run;

/* Run Freqs */
    * These freqs look at the number of genes that a model had the lowest BIC.;
    proc freq data=merged;
        table model /out=freqs;
        run; 

/* Make Gene list of genes with a model better than baseline */
    data nobase;
        set merged;
        where model ^? 'Baseline';
        drop BIC;
        run;

    proc sort data=nobase nodupkey;
        by gene;
        run; *754 genes;

/* Merge on gene symbol */
    proc sort data=DMEL.symbol2fbgn;
        by primary_fbgn;
        run;

    data symmerge;
        retain primary_fbgn symbol model;
        merge nobase (in=in1 rename=(gene=primary_fbgn)) DMEL.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1;
        run;

    proc sort data=symmerge;
        by symbol;
        run;

    data SEM.cegsV_ag_w_BIC_cutoff;
        set symmerge;
        run;

/* Merge to full gene list */
data baseline;
    set SEM.cegsV_gene_list;
    if symbol eq 'Spf45' then flag_baseline = 1;
    else if symbol eq 'Sxl' then flag_baseline = 1;
    else if symbol eq 'Yp2' then flag_baseline = 1;
    else if symbol eq 'dsx' then flag_baseline = 1;
    else if symbol eq 'fl_2_d' then flag_baseline = 1;
    else if symbol eq 'fru' then flag_baseline = 1;
    else if symbol eq 'her' then flag_baseline = 1;
    else if symbol eq 'ix' then flag_baseline = 1;
    else if symbol eq 'snf' then flag_baseline = 1;
    else if symbol eq 'tra' then flag_baseline = 1;
    else if symbol eq 'tra2' then flag_baseline = 1;
    else if symbol eq 'vir' then flag_baseline = 1;
    else flag_baseline = 0;
    run;

proc sort data=baseline;
    by symbol;
    run;

proc sort data=SEM.cegsV_gene2fbgn;
    by gene;
    run;

data fbgn;
    merge baseline (in=in1 rename=(symbol=gene)) SEM.cegsV_gene2fbgn (in=in2);
    by gene;
    if in1;
    keep primary_fbgn flag_baseline;
    run;

proc sort data=SEM.cegsV_ag_w_BIC_cutoff;
    by primary_fbgn;
    run;

proc sort data=fbgn;
    by primary_fbgn;
    run;

data all;
    retain primary_fbgn flag_baseline flag_added;
    merge fbgn (in=in1) SEM.cegsV_ag_w_BIC_cutoff (in=in2);
    by primary_fbgn;
    if in2 then flag_added = 1;
    else flag_added = 0;
    keep primary_fbgn flag_baseline flag_added;
    run;

proc freq data=all;
    tables flag_added;
    tables flag_baseline;
    tables flag_baseline*flag_added;
    run;

data SEM.cegsV_ag_w_flags_bic12;
    set all;
    run;

/* Clean Up */
proc datasets;
    delete merged;
    delete models;
    delete mins;
    delete mins2;
    delete nobase;
    delete freqs;
    delete freqs2;
    delete freqPath;
    delete all;
    delete baseline;
    delete fbgn;
    run; quit;
