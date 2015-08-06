/*******************************************************************************
* Filename: Makefile_dspr_sem_adding_genes.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Summarize adding genes models results.
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* SEM Adding Genome-Wide Gene Level Gene Covariance */
    /* Combine BIC from adding all genes */
        * For each gene, add the baseline model BIC. Format dataset
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/dspr_adding_genes_yp1/sas_data';
        * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/dspr_adding_genes_yp2/sas_data';
        * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/dspr_adding_genes_yp3/sas_data';
        * 
        * INPUT: addgen.*
        *        SEM.dsrp_sex_det_gene_cov_model_bic
        *
        * DATASET: SEM.dsrp_ag_yp1_stack_bic
        *          SEM.dsrp_ag_yp1_sbs_bic
        *          SEM.dsrp_ag_yp2_stack_bic
        *          SEM.dsrp_ag_yp2_sbs_bic
        *          SEM.dsrp_ag_yp3_stack_bic
        *          SEM.dsrp_ag_yp3_sbs_bic
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_adding_genes_yp1.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_adding_genes_yp2.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_adding_genes_yp3.sas';

    /* Identify best model for each gene and see what stands out */
        * 
        * INPUT: SEM.dsrp_ag_yp1_stack_bic
        *        SEM.dsrp_ag_yp2_stack_bic
        *        SEM.dsrp_ag_yp3_stack_bic
        * 
        * BASELINE is always best
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_identify_best_model_yp1.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_identify_best_model_yp2.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_identify_best_model_yp3.sas';
