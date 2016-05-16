/*******************************************************************************
* Filename: tier_import_misspecification.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Luis simulated datasets with differing levels of bias. I ran the
* Bayesian machine will different amounts of misspecification. Import the
* results.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

%macro runLine(line);
    /* Import the Bayesian output and clean up the dataset */
        filename infile "!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/&line._misspecification_results_sim_summary.csv";
        proc import datafile=infile out=CEGS.&line._misspecification dbms=csv replace;
            getnames=yes;
            run;
%mend;
%runLine(r361);
%runLine(r332);
%runLine(r365);
%runLine(r101);
