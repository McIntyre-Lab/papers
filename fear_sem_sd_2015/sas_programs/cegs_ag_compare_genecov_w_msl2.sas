/*******************************************************************************
* Filename: cegs_ag_compare_genecov_w_msl2.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Are the adding gene results different when removing dsx from the
* moddel.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Pull out added genes and compare */
    * While there are similarities in the the counts above, there is clearly a
    * difference in which models were the most likely. Now I want to check and
    * verify that the same genes are added to the expanded network in both
    * models.
    ;

    data genecovExpanded;
        set SEM.cegsV_ag_w_flags_bic12;
        rename flag_expanded = flag_expanded_gcov;
        keep primary_fbgn symbol flag_expanded;
        run;

    data msl2Expanded;
        set SEM.cegsV_ag_w_flags_bic12_w_msl2;
        rename flag_expanded = flag_expanded_mcov;
        keep primary_fbgn symbol flag_expanded;
        run;

    proc sort data=genecovExpanded;
        by primary_fbgn;
        run;

    proc sort data=msl2Expanded;
        by primary_fbgn;
        run;

    data geneMerge;
        merge genecovExpanded (in=in1) msl2Expanded (in=in2);
        by primary_fbgn;
        run;

    proc freq data=geneMerge;
        tables flag_expanded_gcov*flag_expanded_mcov;
        run;

/* Clean up */
    proc datasets nolist;
        delete genecovExpanded;
        delete msl2Expanded;
        delete geneMerge;
        run; quit;

