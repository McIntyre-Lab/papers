/*******************************************************************************
* Filename: cegsV_combine_adding_links_yp2.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Import BIC score from iteratively adding new paths in the SD
* pathway. Format dataset and sort by BIC.
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

libname addgene '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/cegs_adding_links_yp2_fullcov/sas_data';

data addlink;
    set addgene.none;
    drop gene;
    run;

/* merge on paths */
proc sort data=addlink;
    by model;
    run;

proc sort data=SEM.cegsV_al_yp2_model_design_file;
    by model;
    run;

data SEM.cegsV_al_yp2_fullcov_stack_bic;
    retain model modelnum path BIC;
    merge addlink (in=in1) SEM.cegsV_al_yp2_model_design_file (in=in2);
    by model;
    run;

/* Sort by BIC*/
proc sort data=SEM.cegsV_al_yp2_fullcov_stack_bic;
    by BIC;
    run;

proc export data=SEM.cegsV_al_yp2_fullcov_stack_bic outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/cegsV_al_yp2_fullcov_stack_bic.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean Up */
proc datasets nolist;
    delete addlink;
    run; quit;

