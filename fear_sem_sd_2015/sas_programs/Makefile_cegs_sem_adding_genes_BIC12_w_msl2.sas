/*******************************************************************************
* Filename: Makefile_cegs_sem_adding_genes_BIC12_w_msl2.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: From Simulation study, I need to require a BIC differnece of >12
* in order to achieve a 5% TIER. Here I adjust the adding gene results to
* account for this cutoff.
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* DO STUFF ON HPC WITH PYTHON AND SAS SEE DOCUMENTATION */
    * I have a python script that generates sas code for all of the
    * different models. I then use HPC to run all of the SAS programs and
    * output a sas_data set for each gene.
    *
    * SEE: CEGS_sem_adding_genes.xls for more information
    *
    * DATASET: addgen.*
    ;


/* Combine BIC from Adding Genes (Gene Cov With MSL-2) */
    * For each gene, add the baseline model BIC. Format dataset
    *
    * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegs_adding_genes_yp2_w_msl2/sas_data'
    *
    * INPUT: addgen.*
    *        SEM.cegsV_gene_list
    *
    * DATASET: SEM.cegsV_ag_yp2_stack_bic_w_msl2
    *          SEM.cegsV_ag_yp2_sbs_bic_w_msl2
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_genes_nofilter_yp2_w_msl2.sas';

/* Make Model Design file */
    * Models are created using the same iterative process for each genes.
    * So Model 1 will be the same for all genes. I want to create a design
    * file to relate what location the gene was added to the model number.
    * CEGS will be different from DGPR because I used isoforms in DGPR.
    *
    * DATASET: SEM.cegsV_ag_model_design_file
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_ag_create_model_design_file.sas';


/* Export BIC stack with model information (Gene Cov With msl2) */
    * Export the nofilter Yp2 BIC stack with the model paths.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic_w_msl2
    *        SEM.cegsV_ag_model_design_file
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_BIC_w_model_path_w_msl2.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_export_adding_genes_nofilter_yp2_w_msl2.sas';


/* Identify models with BIC 12 less than baseline (Gene Cov With MSL-2) */
    * After doing some simulations, there is a TIER of 10% until we get with a
    * BIC greater in magnitude than 12.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic_w_msl2
    *
    * DATASET: SEM.cegsV_ag_w_BIC_cutoff_w_msl2     334780 models
    *          SEM.cegsV_ag_w_flags_bic12_w_msl2    8810 genes with 0|1 flags for baseline
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_identify_best_model_adding_genes_nofilter_y2_BIC12_w_msl2.sas';


/* Frequency a model was better than baseline (Gene Cov With MSL-2) */
    * Do some models routinely do better than baseline? Are there models
    * that always do worse than baseline? I expect that models corresponding
    * to misspecified parts of the pathway to have many genes. I expect
    * models that are changing the exogenous structure may do worse.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic_w_msl2
    *
    * FILE:
    * !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_freq_a_model_was_better_BIC12_w_msl2.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_dist_models_better_than_baseline_nofilter_yp2_BIC12_w_msl2.sas';


/* Number of models per gene better than baseline (Gene Cov With MSL-2) */
    * Do most genes have multiple models better than baseline? A gene that
    * improves fit for all models may indicate that this gene is really
    * important to the pathway. However, I would not trust the location
    * information. However, if a gene only fits better in a single model
    * than that may indicate that the position is good.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic_w_msl2
    *
    * DATASET:
    * !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_num_models_per_gene_better_than_baseline_BIC12_w_msl2.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_num_models_per_gene_better_than_baseline_nofilter_yp2_BIC12_w_msl2.sas';


/* Compare GeneCov With MSL2 Results */
    * Are the adding gene results different when removing dsx from the model.
    *
    * INPUT: SEM.cegsV_ag_w_bic_cutoff
    *        SEM.cegsV_ag_w_bic_cutoff_w_msl2
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_compare_genecov_w_msl2.sas';
