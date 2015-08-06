/*******************************************************************************
* Filename: Makefile_cegs_sem_adding_new_links.sas
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

/* Gene Covariance Model */
    /* Import new links for CEGS */
        * Import BIC score from iteratively adding new paths in the SD pathway.
        * Format dataset and sort by BIC.
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/adding_newlinks/cegsV_newlink_yp2/sas_data'
        *
        * INPUT: addgen.*
        *
        * DATASET: SEM.cegsV_al_yp2_model_design_file
        *          SEM.cegsV_al_yp2_stack_bic
        *
        * FILES: !MCLAB/cegs_sem_sd_paper/adding_newlinks/cegsV_al_yp2_stack_bic.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_links_yp2.sas';

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

    /* Run the Best model with full covariance model */
        * Comparing the BIC for the Best model (3 Sxl -> Fru BIC = 483.73586) with
        * the same model but with fullcov.
        *
        * INPUT: SEM.cegsV_al_yp2_stack_bic
        *
        * NOTE: BIC for the fullcov is 414.8469
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_run_fullcov_best_adding_links_yp2.sas';

/* Full Covariance Model */
    /* Import new links for CEGS */
        * Import BIC score from iteratively adding new paths in the SD pathway.
        * Format dataset and sort by BIC.
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/adding_newlinks/cegsV_newlink_yp2_fullcof/sas_data'
        *
        * INPUT: addgen.*
        *
        * DATASET: SEM.cegsV_al_yp2_fullcov_stack_bic
        *
        * FILES: !MCLAB/cegs_sem_sd_paper/adding_newlinks/cegsV_al_yp2_fullcov_stack_bic.csv
        ;
        %include %'!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_links_yp2_fullcov.sas';

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
        %include %'!MCLAB/cegs_sem_sd_paper/sas_programs/dspr_cegsV_identify_best_model_newlink_fullcov.sas';

    /* Run the Best model with full covariance model */
        * Comparing the BIC for the Best model (3 Sxl -> Fru BIC = 483.73586) with
        * the same model but with fullcov.
        *
        * INPUT: SEM.cegsV_al_yp2_stack_bic
        *
        * NOTE: BIC for the fullcov is 414.8469
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_run_fullcov_best_adding_links_yp2.sas';
