/*******************************************************************************
* Filename: cegsV_explore_models_nofilter_yp2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Using the dsxNull dataset as a valdiation. Need list of genes
* that were added downstream of DSX.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Figure out what the baseline BIC was */
    * Baseline bic should be the same for all genes. I am grabing the value and
    * putting into a user variable.
    ;
    proc sort data=SEM.cegsV_ag_yp2_stack_bic;
        by gene BIC;
        run;

    data _null_;
        * since sorted by gene the first 38 obs should be the 38 models for the first gene;
        set SEM.cegsV_ag_yp2_stack_bic (obs = 38); 
        where model eq 'Model_Baseline';
        call symput('BIC', BIC);
        run;

    %put _user_;

/* Flag Model if they were less than the Baseline BIC */
    data flag_add;
        set SEM.cegsV_ag_yp2_stack_bic;
        if model eq 'Model_Baseline' then flag_add = 0;
        else if BIC lt &BIC then flag_add = 1;
        else flag_add = 0;
        run;

/* ID models better than baseline */
    proc sort data=flag_add;
        by gene BIC;
        run;

    data improved;
        set flag_add;
        where flag_add = 1;
        dsx_model_rank + 1;
        by gene;
        if first.gene then dsx_model_rank = 1;
        run;
        
/* Any DSX models */
    * All genes that had a DSX model better than baseline. Should have 635
    * genes with model 3 and 604 genes with model 25.
    ;
    data model3;
        set improved;
        if model eq 'Model_3';
        rename gene = symbol;
        rename dsx_model_rank = dsx_model3_rank;
        keep gene dsx_model_rank;
        run;

    data model25;
        set improved;
        if model eq 'Model_25';
        rename gene = symbol;
        rename dsx_model_rank = dsx_model25_rank;
        keep gene dsx_model_rank;
        run;

    proc sort data=model3 nodupkey;
        by symbol;
        run;

    proc sort data=model25 nodupkey;
        by symbol;
        run;

    data flag_all_dsx;
        retain symbol flag_all_dsx_m3 flag_all_dsx_m25;
        merge SEM.cegsV_gene_list (in=in1) model3 (in=in2) model25 (in=in3);
        by symbol;
        if in2 then flag_all_dsx_m3 = 1;
        else flag_all_dsx_m3 = 0;
        if in3 then flag_all_dsx_m25 = 1;
        else flag_all_dsx_m25 = 0;
        run;

    proc freq data=flag_all_dsx;
        tables flag_all_dsx_m3*flag_all_dsx_m25;
        run; * 570 genes had both model 3 and 25, 65 genes had model 3, and 34 genes had model 25;

/* DSX models were best model */
    * Create a list of genes where the DSX model were the best model for that
    * gene. Should have 31 genes with model 3 and 18 genes with model 25.
    ;

    proc sort data=improved;
        by gene BIC;
        run;

    data best_model;
        set improved;
        by gene;
        if first.gene then best = model;
        drop dsx_model_rank;
        run;

    data best_3 best_25;
        set best_model;
        if best eq 'Model_3' then output best_3;
        else if best eq 'Model_25' then output best_25;
        rename gene = symbol;
        keep gene;
        run;

    proc sort data=best_3;
        by symbol;
        run;

    proc sort data=best_25;
        by symbol;
        run;

    data flag_best_dsx;
        merge SEM.cegsV_gene_list (in=in1) best_3 (in=in2) best_25 (in=in3);
        by symbol;
        if in2 then flag_best_dsx_m3 = 1;
        else flag_best_dsx_m3 = 0;
        if in3 then flag_best_dsx_m25 = 1;
        else flag_best_dsx_m25 = 0;
        run;

    proc freq data=flag_best_dsx;
        tables flag_best_dsx_m3*flag_best_dsx_m25;
        run; * 13 genes with model 3, 18 genes with model 25;

/* Merge Flags and add gene symbol */
    proc sort data=flag_all_dsx;
        by symbol;
        run;

    proc sort data=flag_best_dsx;
        by symbol;
        run;

    data merged;
        merge flag_all_dsx (in=in1) flag_best_dsx (in=in2);
        by symbol;
        rename symbol = primary_fbgn;
        run;

    proc sort data=merged;
        by primary_fbgn;
        run;

    proc sort data=DMEL.symbol2fbgn;
        by primary_fbgn;
        run;

    data SEM.cegsV_ag_yp2_flag_ds_dsx;
        retain symbol primary_fbgn;
        merge merged (in=in1) DMEL.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1 and not in2 then symbol = primary_fbgn;
        else if in1;
        run;

/* Export dataset */
    proc export data=SEM.cegsV_ag_yp2_flag_ds_dsx 
        outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_flag_ds_dsx.csv'
        dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
proc datasets nolist;
    delete all_dsx;
    delete best_25;
    delete best_3;
    delete best_model;
    delete flag_add;
    delete flag_all_dsx;
    delete flag_best_dsx;
    delete freq;
    delete improved;
    delete merged;
    delete model25;
    delete model3;
    run; quit;
