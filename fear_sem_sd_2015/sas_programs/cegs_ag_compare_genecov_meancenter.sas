/*******************************************************************************
* Filename: cegs_ag_compare_genecov_meancenter.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Are the adding gene results different when using the gene
* covariance and full covariance models.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create local instance of gene and full covariance results */
    data genecov;
        set SEM.cegsV_ag_w_bic_cutoff;
        rename flag_expanded = flag_expanded_gcov;
        rename flag_most_likely = flag_most_likely_gcov;
        keep primary_fbgn symbol model flag_expanded flag_most_likely;
        run;

    data mc;
        set SEM.cegsV_ag_w_bic_cutoff_mc;
        rename flag_expanded = flag_expanded_mccov;
        rename flag_most_likely = flag_most_likely_mccov;
        keep primary_fbgn symbol model flag_expanded flag_most_likely;
        run;

/* Merge two data sets together */
    proc sort data=genecov;
        by primary_fbgn model;
        run;

    proc sort data=mc;
        by primary_fbgn model;
        run;

    data merged;
        merge genecov (in=in1) mc (in=in2);
        by primary_fbgn model;
        run;

/* Compare results using freqs */
    proc freq data=merged;
        tables flag_expanded_gcov*flag_expanded_mccov /out=expanded;
        tables flag_most_likely_gcov*flag_most_likely_mccov /out=mostlikely;
        run;

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

    data mcExpanded;
        set SEM.cegsV_ag_w_flags_bic12_mc;
        rename flag_expanded = flag_expanded_mccov;
        keep primary_fbgn symbol flag_expanded;
        run;

    proc sort data=genecovExpanded;
        by primary_fbgn;
        run;

    proc sort data=mcExpanded;
        by primary_fbgn;
        run;

    data geneMerge;
        merge genecovExpanded (in=in1) mcExpanded (in=in2);
        by primary_fbgn;
        run;

    proc freq data=geneMerge;
        tables flag_expanded_gcov*flag_expanded_mccov;
        run;

