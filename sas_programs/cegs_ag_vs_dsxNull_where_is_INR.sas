/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/


data ag;
    set SEM.cegsV_ag_nofilter_yp2_stack_bic;
    where gene = 'FBgn0013984';
    run;

/* Model 1 and Model 23 are better than baseline */

