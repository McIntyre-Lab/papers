/*******************************************************************************
* Filename: cegsV_ag_sex_det_explore_dsx.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I am focusing on DSX for the paper, so I need to look at where
* it fits in best.
*
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

data dsx;
    set SEM.cegsV_ag_sex_det_genes_stack_bic;
    where gene eq 'dsx';
    run;

proc sort data=dsx;
    by BIC;
    run;
