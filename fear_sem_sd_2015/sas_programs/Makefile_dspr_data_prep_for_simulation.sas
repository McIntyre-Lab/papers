/*******************************************************************************
* Filename: Makefile_dspr_data_prep_for_simulation.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Prepare the data for doing simulations 
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Create null list of genes */
    * Create a list of genes that are not associated with the sex hierarchy.
    * Merge on the added genes and keep only those genes not added to the SD
    * pathway.
    *
    * INPUT: SEM.dsrp_stack_gene_level_sym
    * 
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/dspr_all_distributions.pdf
    *       !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/dspr_all_mean_and_variances.csv    n=7,423
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dspr_simulation_make_all_list.sas';
