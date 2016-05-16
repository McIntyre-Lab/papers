/*******************************************************************************
* Filename: tier_import_r101_misspecification.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Luis simulated datasets with differing levels of bias. I ran the Bayesian
* machine will different amounts of misspecification. Import the results.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Import the Bayesian output and clean up the dataset */
    proc import datafile='!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/output/old/r101_simulateddata_forTIER_FigureNOAIPoisson_q_true_results_summary.csv' out=CEGS.r101_fb557_miss dbms=csv replace;
        getnames=yes;
        run;

