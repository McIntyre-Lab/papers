/*******************************************************************************
* Filename: dspr_cegsV_identify_best_model_newlink.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: In order to refine which relationships to try to add to the
* network, I am only taking putative paths that are present in both DSPR and
* CEGS dataset. I am also arbitraily requiring that the BIC score be less than
* BASELINE-2.
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Grab BIC for Baselines */
data cbaseline;
    set SEM.cegsV_al_yp2_fullcov_stack_bic;
    where modelnum = 0;
    run;

data dbaseline;
    set SEM.dspr_al_yp2_fullcov_stack_bic;
    where modelnum = 0;
    run;

data _null_;
    set cbaseline;
    call symput('cbase', BIC);
    run;

data _null_;
    set dbaseline;
    call symput('dbase', BIC);
    run;

/* Pull added paths that have BIC with less than baseine BIC - 2 */
data cegs;
    set SEM.cegsV_al_yp2_fullcov_stack_bic;
    where BIC le &cbase - 2 ;
    rename BIC = cBIC;
    run; * 29 paths;

data dspr;
    set SEM.dspr_al_yp2_fullcov_stack_bic;
    where BIC le &dbase - 2 ;
    rename BIC = dBIC;
    run; * 51 paths;

/* Merge by Path, keep paths that are in both datasets */
proc sort data=cegs;
    by path;
    run;

proc sort data=dspr;
    by path;
    run;

data merged;
    merge cegs (in=in1) dspr (in=in2);
    by path;
    if in1 and in2;
    drop model modelnum;
    run; * 16 paths;

/* output list for michelle */
proc sort data=merged;
    by cBIC;
    run;

proc export data=merged outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/putative_paths_fullcov.csv' dbms=csv replace;
    putnames=yes;
    run;
