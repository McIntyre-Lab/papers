/*******************************************************************************
* Filename: cegsV_identify_best_model_adding_genes_nofilter_y2_BIC12_line_subset.sas
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
        set SEM.cegsV_ag_yp2_stack_bic_subset ;
        where model eq 'Model_Baseline';
        if _n_ eq 1 then do;
        call symput('Baseline', BIC);
        end;
        run;

    %put _user_;

/* How many models*gene imporved fit */
    data models;
        set SEM.cegsV_ag_yp2_stack_bic_subset ;
        if BIC le &BASELINE - 12 then flag_expanded = 1;
        else flag_expanded = 0;
        run; 
        
    proc freq data=models;
        table flag_expanded;
        run; * 9,548 models better than baseline -- 325,323 worse than baseline;

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
        run; *770 genes;

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

    data SEM.cegsV_ag_w_BIC_cutoff_subset;
        set symmerge;
        keep primary_fbgn symbol model gene BIC flag_expanded flag_most_likely;
        run;

/* Reduce Models to Gene list */
    proc sort data=SEM.cegsV_ag_w_BIC_cutoff_subset;
        by primary_fbgn BIC;
        run;

    data gene;
        set SEM.cegsV_ag_w_BIC_cutoff_subset;
        by primary_fbgn;
        if first.primary_fbgn then output;
        drop model BIC flag_most_likely;
        run;

    data SEM.cegsV_ag_w_flags_bic12_subset;
        set gene;
        run;

/* Clean Up */
proc datasets nolist;
    delete models;
    delete mins;
    delete merged;
    delete gene;
    delete symmerge;
    run; quit;
