/*******************************************************************************
* Filename: Makefile_cegs_sem_adding_genes_BIC12_genecov.sas
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


/* Combine BIC from Adding Genes (Gene Cov) */
    /* Yp1 */
        * For each gene, add the baseline model BIC. Format dataset
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegs_adding_genes_yp1/sas_data'
        *
        * INPUT: addgen.*
        *        SEM.cegsV_gene_list
        *
        * DATASET: SEM.cegsV_ag_yp1_stack_bic
        *          SEM.cegsV_ag_yp1_sbs_bic
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_genes_nofilter_yp1.sas';

    /* Yp2 */
        * For each gene, add the baseline model BIC. Format dataset
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegs_adding_genes_yp2/sas_data'
        *
        * INPUT: addgen.*
        *        SEM.cegsV_gene_list
        *
        * DATASET: SEM.cegsV_ag_yp2_stack_bic
        *          SEM.cegsV_ag_yp2_sbs_bic
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_genes_nofilter_yp2.sas';

    /* Yp3 */
        * For each gene, add the baseline model BIC. Format dataset
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegs_adding_genes_yp3/sas_data'
        *
        * INPUT: addgen.*
        *        SEM.cegsV_gene_list
        *
        * DATASET: SEM.cegsV_ag_yp3_stack_bic
        *          SEM.cegsV_ag_yp3_sbs_bic
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_combine_adding_genes_nofilter_yp3.sas';


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


/* Export BIC stack with model information (Gene Cov) */
    * Export the nofilter Yp2 BIC stack with the model paths.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic
    *        SEM.cegsV_ag_model_design_file
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_BIC_w_model_path.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_export_adding_genes_nofilter_yp2.sas';


/* Identify models with BIC 12 less than baseline (Gene Cov) */
    * After doing some simulations, there is a TIER of 10% until we get with a
    * BIC greater in magnitude than 12.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic
    *
    * DATASET: SEM.cegsV_ag_w_BIC_cutoff     334780 models
    *          SEM.cegsV_ag_w_flags_bic12    8810 genes with 0|1 flags for baseline
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_identify_best_model_adding_genes_nofilter_y2_BIC12_v2.sas';


/* Frequency a model was better than baseline (Gene Cov) */
    * Do some models routinely do better than baseline? Are there models
    * that always do worse than baseline? I expect that models corresponding
    * to misspecified parts of the pathway to have many genes. I expect
    * models that are changing the exogenous structure may do worse.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic
    *
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_freq_a_model_was_better_BIC12.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_dist_models_better_than_baseline_nofilter_yp2_BIC12.sas';


/* Number of models per gene better than baseline (Gene Cov) */
    * Do most genes have multiple models better than baseline? A gene that
    * improves fit for all models may indicate that this gene is really
    * important to the pathway. However, I would not trust the location
    * information. However, if a gene only fits better in a single model
    * than that may indicate that the position is good.
    *
    * INPUT: SEM.cegsV_ag_yp2_stack_bic
    *
    * DATASET: !MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/cegsV_ag_yp2_num_models_per_gene_better_than_baseline_BIC12.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_num_models_per_gene_better_than_baseline_nofilter_yp2_BIC12.sas';


/* Compare GeneCov FullCov results */
    * Are the adding gene results different when using the gene covariance and
    * full covariance models.
    *
    *       ### COMPARISON OF MODELS ADDED ###
    *       flag_expanded_gcov flag_expanded_fcov
    *
    *       Frequency|
    *       Percent  |
    *       Row Pct  |
    *       Col Pct  |       0|       1|  Total
    *       ---------+--------+--------+
    *              0 | 321716 |    499 | 322215
    *                |  96.10 |   0.15 |  96.25
    *                |  99.85 |   0.15 |
    *                |  99.76 |   4.06 |
    *       ---------+--------+--------+
    *              1 |    775 |  11790 |  12565
    *                |   0.23 |   3.52 |   3.75
    *                |   6.17 |  93.83 |
    *                |   0.24 |  95.94 |
    *       ---------+--------+--------+
    *       Total      322491    12289   334780
    *                   96.33     3.67   100.00
    *
    *
    *       ### COMPARISON OF MOST LIKELY MODELS ADDED ###
    *       flag_most_likely_gcov flag_most_likely_fcov
    *
    *       Frequency|
    *       Percent  |
    *       Row Pct  |
    *       Col Pct  |       0|       1|  Total
    *       ---------+--------+--------+
    *              0 | 333473 |    553 | 334026
    *                |  99.61 |   0.17 |  99.77
    *                |  99.83 |   0.17 |
    *                |  99.89 |  59.78 |
    *       ---------+--------+--------+
    *              1 |    382 |    372 |    754
    *                |   0.11 |   0.11 |   0.23
    *                |  50.66 |  49.34 |
    *                |   0.11 |  40.22 |
    *       ---------+--------+--------+
    *       Total      333855      925   334780
    *                   99.72     0.28   100.00
    *
    *
    *       ### COMPARISON OF GENES ADDED ###
    *       flag_expanded_gcov flag_expanded_fcov
    *
    *       Frequency|
    *       Percent  |
    *       Row Pct  |
    *       Col Pct  |       0|       1|  Total
    *       ---------+--------+--------+
    *              0 |   7885 |    171 |   8056
    *                |  89.50 |   1.94 |  91.44
    *                |  97.88 |   2.12 |
    *                | 100.00 |  18.49 |
    *       ---------+--------+--------+
    *              1 |      0 |    754 |    754
    *                |   0.00 |   8.56 |   8.56
    *                |   0.00 | 100.00 |
    *                |   0.00 |  81.51 |
    *       ---------+--------+--------+
    *       Total        7885      925     8810
    *                   89.50    10.50   100.00
    *
    *
    * INPUT: SEM.cegsV_ag_w_bic_cutoff
    *        SEM.cegsV_ag_fullcov_w_BIC_cutoff
    *
    * DATASET:
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_compare_genecov_fullcov.sas';

/* Compare GeneCov LineSubset results */
    * Are the adding gene results different when using a subset of the lines.
    *
    * INPUT: SEM.cegsV_ag_w_bic_cutoff
    *        SEM.cegsV_ag_w_bic_cutoff_subset
    *
    *       flag_expanded_gcov flag_expanded_scov
    *
    *       Frequency|
    *       Percent  |
    *       Row Pct  |
    *       Col Pct  |       0|       1|  Total
    *       ---------+--------+--------+
    *              0 |   7762 |    294 |   8056
    *                |  88.10 |   3.34 |  91.44
    *                |  96.35 |   3.65 |
    *                |  96.54 |  38.18 |
    *       ---------+--------+--------+
    *              1 |    278 |    476 |    754
    *                |   3.16 |   5.40 |   8.56
    *                |  36.87 |  63.13 |
    *                |   3.46 |  61.82 |
    *       ---------+--------+--------+
    *       Total        8040      770     8810
    *                   91.26     8.74   100.00
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_compare_genecov_subset.sas';

/* Compare GeneCov No DSX Results */
    * Are the adding gene results different when removing dsx from the model.
    *
    * INPUT: SEM.cegsV_ag_w_bic_cutoff
    *        SEM.cegsV_ag_w_bic_cutoff_no_dsx
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegs_ag_compare_genecov_no_dsx.sas';


/* Export genes for enrichment analysis */
    proc sort data=SEM.cegsV_ag_w_bic_cutoff;
        by symbol bic;
        run;

    proc export data=SEM.cegsV_ag_w_bic_cutoff outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/adding_genes/gene_list_bic_cutoff.csv' dbms=csv replace;
        putnames=yes;
        run;
