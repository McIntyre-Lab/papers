/*******************************************************************************
* Filename: cegsV_explore_models_nofilter_yp2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Do some models routinely do better than baseline? Are there
* models that always do worse than baseline? I expect that models corresonding
* to mispecified parts of the pathway to have many genes. I expect models that
* are changing the exogenous structure may do worse. 
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
        * since I sorted by gene the first 38 obs should be all 38 models for the first gene;
        set SEM.cegsV_ag_yp2_stack_bic (obs = 38); 
        where model eq 'Model_Baseline';
        call symput('BIC', BIC);
        run;

    %put _user_;

/* Flag Model if they were better than baseline */
    data flag_add;
        set SEM.cegsV_ag_yp2_stack_bic ;
        if model eq 'Model_Baseline' then flag_add = 0;
        else if BIC lt &BIC then flag_add = 1;
        else flag_add = 0;
        run;

/* Distribution of models better than baseline */
    data dist;
        set flag_add;
        where flag_add = 1;
        run;

    proc freq data=dist;
        tables model /out=freqdist;
        run;

    * Merge on path information;
    proc sort data=freqdist;
        by model;
        run;
    
    proc sort data=SEM.cegsV_ag_model_design_file;
        by model;
        run;

    data freqdist_w_path;
        merge freqdist (in=in1) SEM.cegsV_ag_model_design_file (in=in2);
        by model;
        drop model percent;
        run;

    data freqdist_w_path2;
        retain modelnum path;
        set freqdist_w_path;
        if modelnum eq 0 then delete;
        rename modelnum = model;
        label count = ' ';
        run;

    proc sort data=freqdist_w_path2;
        by model;
        run;

/* Export dataset */
    proc export data=freqdist_w_path2 
        outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_freq_a_model_was_better.csv' 
        dbms=csv replace;
        putnames=yes;
        run;

/* Clean up */
proc datasets nolist;
    delete dist;
    delete flag_add;
    delete freqdist;
    delete freqdist_w_path;
    delete freqdist_w_path2;
    run; quit;

