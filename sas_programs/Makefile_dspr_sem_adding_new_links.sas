/*******************************************************************************
* Filename: Makefile_dspr_sem_adding_new_links.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: The goal of this part is to take an exisiting network structure
* and add new relationships (paths) between the genes already present in the
* network.
*******************************************************************************/

/* Libraries */

libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Gene Covariance model */
    /* Import new links for DSPR */
        * Import BIC score from iteratively adding new paths in the SD pathway.
        * Format dataset and sort by BIC.
        *
        * libname addgene '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/dspr_adding_links_yp2/sas_data';
        *
        * INPUT: addgen.*
        *
        * DATASET: SEM.dspr_al_yp2_model_design_file
        *          SEM.dspr_al_yp2_stack_bic
        *
        * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/dspr_al_yp2_stack_bic.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dspr_combine_adding_links_yp2.sas';

    /* Combine DSPR and CEGS Models to identify putative new relationships */
        * In order to refine which relationships to try to add to the network,
        * I am only taking putative paths that are present in both DSPR and
        * CEGS dataset. I am also arbitraily requiring that the BIC score be
        * less than BASELINE-2.
        *
        * INPUT: SEM.cegsV_al_yp2_stack_bic
        *        SEM.dspr_al_yp2_stack_bic
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/putative_paths.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dspr_cegsV_identify_best_model_newlink.sas';

/* Full Covariance model */
    /* Import new links for DSPR */
        * Import BIC score from iteratively adding new paths in the SD pathway.
        * Format dataset and sort by BIC.
        *
        * libname addgene '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/dspr_adding_links_yp2_fullcov/sas_data';
        *
        * INPUT: addgen.*
        *
        * DATASET: SEM.dspr_al_yp2_fullcov_stack_bic
        *
        * FILES: !MCLAB/cegs_sem_sd_paper/adding_newlinks/dspr_al_yp2_stack_bic.csv
        ;
        %include %'!MCLAB/cegs_sem_sd_paper/sas_programs/dspr_combine_adding_links_yp2_fullcov.sas';

    /* Combine DSPR and CEGS Models to identify putative new relationships */
        * In order to refine which relationships to try to add to the network,
        * I am only taking putative paths that are present in both DSPR and
        * CEGS dataset. I am also arbitraily requiring that the BIC score be
        * less than BASELINE-2.
        *
        * INPUT: SEM.cegsV_al_yp2_fullcov_stack_bic
        *        SEM.dspr_al_yp2_fullcov_stack_bic
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_newlinks/putative_paths_fullcov.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dspr_cegsV_identify_best_model_newlink.sas';

