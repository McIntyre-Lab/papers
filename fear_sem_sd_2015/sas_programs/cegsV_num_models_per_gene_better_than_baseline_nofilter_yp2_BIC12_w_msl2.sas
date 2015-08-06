/*******************************************************************************
* Filename:
* cegsV_num_models_per_gene_better_than_baseline_nofilter_yp2_BIC12_w_msl2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Do most genes have multiple models better than baseline? A gene
* that improves fit for all models may indicate that this gene is really
* important to the pathway. However, I would not trust the loction information.
* However, if a gene only fits better in a single model than that may indicate
* that the position is good.
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

    proc sort data=SEM.cegsV_ag_yp2_stack_bic_w_msl2;
        by gene BIC;
        run;

    data _null_;
        * since sorted by gene the first 38 obs should be the 38 models for the first gene;
        set SEM.cegsV_ag_yp2_stack_bic_w_msl2 (obs = 38); 
        where model eq 'Model_Baseline';
        call symput('BIC', BIC);
        run;

    %put _user_;

/* Flag Model if they were less than the Baseline BIC */
    data flag_add;
        set SEM.cegsV_ag_yp2_stack_bic_w_msl2;
        if model eq 'Model_Baseline' then flag_add = 0;
        else do;
        diff = &BIC - BIC;
        if diff ge 12 then flag_add = 1;
        else flag_add = 0;
        end;
        drop diff;
        run;

/* Distribution of the number of models better by gene */
    * If a gene is added, is only a single model flagged, or is it typically
    * multiple models. 
    ;

    proc sort data=flag_add;
        by gene;
        run;

    proc means data=flag_add noprint;
        by gene;
        output out=sums sum(flag_add)=sum_flag_add;
        run;

    proc freq data=sums;
        table sum_flag_add /out=freq;
        run;
        
/* Export dataset */
proc export data=freq 
    outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_num_models_per_gene_better_than_baseline_BIC12_w_msl2.csv'
    dbms=csv replace;
    putnames=yes;
    run;

/* Clean up */
proc datasets nolist;
    delete flag_add;
    delete freq;
    delete sums;
    run; quit;
