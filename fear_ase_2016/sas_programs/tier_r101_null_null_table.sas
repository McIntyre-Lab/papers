/*******************************************************************************
* Filename: tier_r101_null_null_table.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Like the Luis paper, I want to create a table with the TIER
* percentages for the NoBias NoAI simulated data.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

proc freq data=noai_nobias;
    tables flag_AI_binomial_test flag_AI_NB_p_random_DNA
    flag_AI_PG_q_random_DNA flag_AI_PG_q4 flag_AI_PG_q5 flag_AI_PG_q6 ;
    run;

data tmp;
    set noai_nobias;
    if flag_AI_PG_q4 = 1 and flag_AI_PG_q5 = 1 and flag_AI_PG_q6 = 1 then flag_all = 1;
    else flag_all = 0;
    run;

proc freq data=tmp;
    tables flag_all;
    run;

proc datasets nolist;
    delete tmp;
    run; quit;

