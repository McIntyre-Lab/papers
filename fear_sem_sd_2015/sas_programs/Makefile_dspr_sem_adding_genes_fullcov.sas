/*******************************************************************************
* Filename: Makefile_dspr_sem_adding_genes_fullcov.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Run adding genes models using the full covarance baseline model.
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* SEM Adding Genome-Wide Gene Level Full Covariance */
    /* Combine BIC from adding all genes */
        * For each gene, add the baseline model BIC. Format dataset
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/dspr_adding_genes_yp2_fullcov/sasdata';
        * 
        * INPUT: addgen.*
        *        SEM.dsrp_gene_list
        *
        * DATASET: 
        *          SEM.dsrp_ag_yp2_fullcov_stack_bic
        *          SEM.dsrp_ag_yp2_fullcov_sbs_bic
        ;
        %include %'!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_adding_genes_yp2_fullcov.sas';

    /* Identify best model for each gene and see what stands out */
        * 
        * INPUT: SEM.dsrp_ag_yp2_fullcov_stack_bic
        * 
        * BASELINE is always best
        ;
        %include %'!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_identify_best_model_yp2_fullcov.sas';
