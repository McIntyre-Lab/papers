/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

proc means data=SEM.dsrp_ag_yp2_fullcov_stack_bic noprint;
    by gene;
    var BIC;
    output out=mins min=min;
    run;

data mins2;
    set mins;
    rename min=BIC;
    keep gene min;
    run;

proc sort data=mins2;
    by gene BIC;
    run;

proc sort data=SEM.dsrp_ag_yp2_fullcov_stack_bic;
    by gene BIC;
    run;

data merged;
    merge SEM.dsrp_ag_yp2_fullcov_stack_bic  (in=in1) mins2 (in=in2);
    by gene BIC;
    if in1 and in2;
    run;

proc sort data=merged;
    by model;
    run;

proc freq data=merged;
    table model;
    run;
