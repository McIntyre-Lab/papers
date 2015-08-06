/*******************************************************************************
* Filename: cegsV_identify_best_model_adding_genes_nofilter_y2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: For each gene figure out the lowest BIC score and create a list
* of genes that improved overall model fit above baseline.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

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
        run;

    proc sort data=merged;
        by model;
        run;

/* Run Freqs */
    * These freqs look at the number of genes that a model had the lowest BIC.;
    proc freq data=merged;
        table model /out=freqs;
        run;

    * Merge on path information;
    proc sort data=freqs;
        by model;
        run;

    proc sort data=SEM.cegsV_ag_model_design_file;
        by model;
        run;

    data freqPath;
        retain model path;
        merge freqs (in=in1) SEM.cegsV_ag_model_design_file (in=in2);
        by model;
        if modelnum = '.' then modelnum=0;
        label count = ' ';
        drop percent;
        run;

    * sort and export;
    proc sort data=freqPath;
        by modelnum;
        run;

    data freqPath;
        set freqPath;
        drop modelnum;
        run;

    proc export data=freqPath 
            outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_freq_a_model_was_best.csv'
            dbms=csv replace;
        putnames=yes;
        run;

/* Make Gene list of genes with a model better than baseline */
    data nobase;
        set merged;
        where model ^? 'Baseline';
        drop BIC;
        run;

    proc sort data=nobase nodupkey;
        by gene;
        run; *1390 genes, matches sum from freqs;

    data SEM.cegsV_ag_yp2_added_genes;
        set nobase;
        run;

/* Clean Up */
proc datasets;
    delete merged;
    delete mins;
    delete mins2;
    delete nobase;
    delete freqs;
    delete freqs2;
    delete freqPath;
    run; quit;
