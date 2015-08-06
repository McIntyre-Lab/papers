/*******************************************************************************
* Filename: cegsV_explore_models_nofilter_yp2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: It would be useful to know what genes are downstream of fru.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sas_data';
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
    data improved;
        set flag_add;
        where flag_add = 1;
        run;
        
/* Any fru models */
    * All genes that had a Fru model (Model_4) better than baseline. 
    ;
    data model4 ;
        set improved;
        if model eq 'Model_4' then output model4;
        rename gene = symbol;
        keep gene;
        run;

    proc sort data=model4 nodupkey;
        by symbol;
        run;

    data flag_all_fru;
        merge SEM.cegsV_gene_list (in=in1) model4 (in=in2);
        by symbol;
        if in2 then flag_all_fru = 1;
        else flag_all_fru = 0;
        run;

    proc freq data=flag_all_fru;
        tables flag_all_fru;
        run; * 804 genes had model 4;

/* Fru models were best model */
    * Create a list of genes where the Fru model were the best model for that
    * gene. 
    ;
    proc sort data=improved;
        by gene BIC;
        run;

    data best_model;
        set improved;
        by gene;
        if first.gene then best = model;
        run;

    data best_4;
        set best_model;
        if best eq 'Model_4' then output best_4;
        rename gene = symbol;
        keep gene;
        run;

    proc sort data=best_4;
        by symbol;
        run;

    data flag_best_fru;
        merge SEM.cegsV_gene_list (in=in1) best_4 (in=in2);
        by symbol;
        if in2 then flag_best_fru = 1;
        else flag_best_fru = 0;
        run;

    proc freq data=flag_best_fru;
        tables flag_best_fru;
        run; * 9 genes had model 4 as best model;

/* Merge Flags and add gene symbol */
    proc sort data=flag_all_fru;
        by symbol;
        run;

    proc sort data=flag_best_fru;
        by symbol;
        run;

    data merged;
        merge flag_all_fru (in=in1) flag_best_fru (in=in2);
        by symbol;
        rename symbol = primary_fbgn;
        run;

    proc sort data=merged;
        by primary_fbgn;
        run;

    proc sort data=DMEL.symbol2fbgn;
        by primary_fbgn;
        run;

    data SEM.cegsV_ag_yp2_flag_ds_fru;
        retain symbol primary_fbgn;
        merge merged (in=in1) DMEL.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1 and not in2 then symbol = primary_fbgn;
        else if in1;
        run;

/* Export dataset */
    proc export data=SEM.cegsV_ag_yp2_flag_ds_fru 
        outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_flag_ds_fru.csv'
        dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
proc datasets nolist;
    delete all_fru;
    delete best_4;
    delete best_model;
    delete flag_add;
    delete flag_all_fru;
    delete flag_best_fru;
    delete freq;
    delete improved;
    delete merged;
    delete model4;
    run; quit;
