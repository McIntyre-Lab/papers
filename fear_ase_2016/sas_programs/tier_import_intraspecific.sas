/*******************************************************************************
* Filename: tier_import_intraspecific.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Using intraspecific data simulated from r101, I ran the Bayesian
* analysis.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

proc import datafile='!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/output/r101_simulated_NOAI_NOBIASRequal1Poisson_results_summary.csv' out=noai_nobias dbms=csv replace;
    getnames=yes;
    guessingrows=99999;
    run;


/*
proc import datafile='!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/output/simulated_NOAI_YESBIASRequal1Poisson_results_summary.csv' out=noai_bias dbms=csv replace;
    getnames=yes;
    run;


proc import datafile='!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/output/simulated_YESAI_YESBIASRequal1d5Poisson_results_summary.csv' out=ai_bias dbms=csv replace;
    getnames=yes;
    run;
*/
