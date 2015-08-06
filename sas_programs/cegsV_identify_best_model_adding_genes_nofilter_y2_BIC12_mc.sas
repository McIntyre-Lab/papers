/*******************************************************************************
* Filename: cegsV_identify_best_model_adding_genes_nofilter_y2_BIC12_mc.sas
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
        set SEM.cegsV_ag_yp2_stack_bic_mc ;
        where model eq 'Model_Baseline';
        if _n_ eq 1 then do;
        call symput('Baseline', BIC);
        end;
        run;

/* How many models*gene imporved fit */
    data models;
        set SEM.cegsV_ag_yp2_stack_bic_mc ;
        if BIC le &BASELINE - 12 then flag_expanded = 1;
        else flag_expanded = 0;
        run; 
        
    proc freq data=models;
        table flag_expanded;
        run; * 12,565 models better than baseline;

/* For each gene figure out the lowest BIC */
    proc sort data=models;
        by gene bic;
        run;

    data mins;
        format gene $300.;
        set models;
        by gene;
        if first.gene and flag_expanded = 1 then flag_most_likely = 1;
        else flag_most_likely = 0;
        run;

/* Run Freqs */
    * These freqs look at the number of genes that a model had the lowest BIC.;
    proc freq data=mins(where=(flag_most_likely = 1));
        table model;
        run; 

/* Make Gene list of genes with a model better than baseline */
    proc freq data=mins;
        table flag_most_likely;
        run; *754 genes;

/* Merge on primary_fbgn */
    proc sort data=mins;
        by gene;
        run;

    proc sort data=SEM.cegsV_gene2fbgn;
        by gene;
        run;
    
    data merged;
        retain primary_fbgn;
        merge mins (in=in1) SEM.cegsV_gene2fbgn (in=in2);
        by gene;
        if in1;
        run;

/* Merge on gene symbol */
    proc sort data=DMEL.symbol2fbgn;
        by primary_fbgn;
        run;

    proc sort data=merged;
        by primary_fbgn;
        run;

    data symmerge;
        retain primary_fbgn symbol model;
        merge merged (in=in1) DMEL.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1;
        run;

    proc sort data=symmerge;
        by gene;
        run;

    proc sort data=SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12;
        by primary_fbgn;
        run;

    proc sort data=SEM.cegsV_ag_yp2_flag_ds_fru_BIC12;
        by primary_fbgn;
        run;

    data merged;
        merge symmerge (in=in1)  SEM.cegsV_ag_yp2_flag_ds_dsx_BIC12 (in=in2 rename=(primary_fbgn = gene)) SEM.cegsV_ag_yp2_flag_ds_fru_BIC12 (in=in3 rename=(primary_fbgn = gene));
        by gene;
        if in1 ;
        drop dsx_model3_rank dsx_model25_rank gene;
        run;

    data SEM.cegsV_ag_w_BIC_cutoff_mc;
        set merged;
        if flag_all_fru = 1 then flag_fru_exp = 1;
        else flag_fru_exp = 0;
        if flag_best_fru = 1 then flag_fru_exp_most_likely = 1;
        else flag_fru_exp_most_likely = 0;
        if flag_all_dsx_m3 = 1 or flag_all_dsx_m25 = 1 then flag_dsx_exp = 1;
        else flag_dsx_exp = 0;
        if flag_best_dsx_m3 = 1 or flag_best_dsx_m25 = 1 then flag_dsx_exp_most_likely = 1;
        else flag_dsx_exp_most_likely = 0;
        drop flag_all_dsx_m3 flag_all_dsx_m25 flag_best_dsx_m3 flag_best_dsx_m25 flag_all_fru flag_best_fru;
        run;

/* Reduce Models to Gene list */
    proc sort data=SEM.cegsV_ag_w_BIC_cutoff_mc;
        by primary_fbgn BIC;
        run;

    data gene;
        set SEM.cegsV_ag_w_BIC_cutoff_mc;
        by primary_fbgn;
        if first.primary_fbgn then output;
        drop model BIC flag_most_likely;
        run;

    data SEM.cegsV_ag_w_flags_bic12_mc;
        set gene;
        run;

    proc freq data=SEM.cegsV_ag_w_flags_bic12_mc;
        table flag_fru_exp;
        run;

/* Clean Up */
proc datasets;
    delete models;
    delete mins;
    delete merged;
    delete gene;
    delete symmerge;
    run; quit;
