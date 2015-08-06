/*******************************************************************************
* Filename: Makefile_cegs_data_prep_for_simulation.sas
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

/* Create all list of genes */
    * Create a list of genes that are not associated with the sex hierarchy.
    * Merge on the added genes and keep only those genes not added to the SD
    * pathway.
    *
    * INPUT: SEM.cegsV_by_gene_stack
    * 
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_all_distributions.pdf
    *       !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_all_mean_and_variances.csv    n=8,882
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegV_simulation_make_all_list.sas';

/* Create null list of genes */
    * Create a list of genes that are not associated with the sex hierarchy.
    * Merge on the added genes and keep only those genes not added to the SD
    * pathway.
    *
    * INPUT: SEM.cegsV_ag_yp2_added_genes
    *        SEM.cegsV_gene_list
    *        SEM.cegsV_by_gene_stack
    * 
    * FILE: 
    *       !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_null_distributions.pdf
    *       !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/cegsV_null_mean_and_variances.csv    n=7,419
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegV_simulation_make_null_list.sas';

/* Create distribution of added genes */
    * I am curious to see if the genes added to the SD pathway have a different
    * distribution of means and variances.
    *
    * INPUT: SEM.cegsV_ag_yp2_added_genes
    *        SEM.cegsV_by_gene_stack
    * 
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/simulation/added_distributions.pdf
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegV_simulation_added_distribution.sas';
