/*******************************************************************************
* Filename: qsim_import_pg_results.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Import the PG model results using q = qsim.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname luis '!MCLAB/CEGS_Bayesian_Analysis_Paper/SAS_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

data CEGS.qsim_bayesian_results;
    infile '!MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/output/ase_dataset_for_bayesian_w_qsim_summary.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
    informat line $5. ;
    informat mating_status $1. ;
    informat fusion_id $20. ;
    informat qsim_mean_theta best32. ;
    informat qsim_q025 best32. ;
    informat qsim_q975 best32. ;
    informat Bayesianpvalue_qsim best32.;
    informat flag_AI_qsim best32.;
    format line $5. ;
    format mating_status $1. ;
    format fusion_id $20. ;
    format qsim_mean_theta best12. ;
    format qsim_q025 best12. ;
    format qsim_q975 best12. ;
    format Bayesianpvalue_qsim best12.;
    format flag_AI_qsim best12.;
    input
             line $
             mating_status $
             fusion_id $
             qsim_mean_theta 
             qsim_q025 
             qsim_q975 
             Bayesianpvalue_qsim
             flag_AI_qsim
    ;
    run;
