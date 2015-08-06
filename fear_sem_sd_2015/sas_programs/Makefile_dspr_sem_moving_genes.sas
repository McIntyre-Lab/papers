libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* SEM Moving genes around one-at-a-time */
    * Repeat a basic SEM model moving some of the genes around. Note I am
    * going to select only a single isoform for each gene to use in this
    * model.
    *
    * INPUT: SEM.dsrp_sex_det_sbs_combine_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_yp1.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_yp3.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_spf45_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_vir.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_fl2d.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_her.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_yp2.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_yp2_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/spf45_add_fru_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_fru_only.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_tra_only.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_tra2a_only.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_yp2_only.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_her_only.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_only.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_fru_to_yp2_only.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_w_spf45_to_fru_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_no_spf45_estimates.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_model.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_yp1.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_yp3.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_all_yps.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_estimates.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_spf45_estimates.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_vir.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_fl2d.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_her.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/remove_yp2.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_yp2_estimates.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/spf45_add_fru_estimates.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_fru_only.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_tra_only.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_tra2a_only.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_yp2_only.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_spf45_to_her_only.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_only.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_fru_to_yp2_only.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_w_spf45_to_fru_estimates.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_no_spf45_estimates.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_estimates.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_yp1.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_yp3.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_all_yp3.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/remove_spf45_estimates.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/remove_vir.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/remove_fl2d.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/remove_her.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/remove_yp2.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/spf45_add_fru_estimates.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/spf45_add_fru_estimates.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_spf45_to_fru_only.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_spf45_to_tra_only.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_spf45_to_yp2_only.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_spf45_to_her_only.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_vir_to_tra2_only.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_fru_to_yp2_only.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/move_vir_to_tra2_w_spf45_to_fru_estimates.log
    *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/move_vir_to_tra2_no_spf45_estimates.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_gene_cov.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_gene_cov_yp1.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_gene_cov_yp3.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_gene_cov_all_yps.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_full_cov.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_remove_spf45.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_remove_vir.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_remove_fl2d.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_remove_her.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_remove_yp2.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_spf45_to_yp2.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_spf45_to_fru.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_spf45_to_fru_only.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_spf45_to_tra_only.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_spf45_to_tra2_only.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_spf45_to_yp2_only.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_spf45_to_her_only.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_vir_to_tra2.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_fru_to_yp2.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_vir_to_tra2_move_spf45_to_fru.sas';
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_move_vir_to_tra2_remove_spf45.sas';


/* SEM Moving genes around Genome-Wide Collapsed Isoforms */
    /* Combine BIC from moving all genes */
        * For each gene, add the baseline model BIC. Format dataset
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/adding_genes/newlink/sas_data'
        * 
        * INPUT: addgen.*
        *        SEM.dsrp_sex_det_gene_cov_model_bic
        *
        * DATASET: SEM.dsrp_add_newlink_stack_bic
        *          SEM.dsrp_add_newlink_sbs_bic
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_combine_newlink.sas';

    /* Identify best model for each gene and see what stands out */
        * 
        * INPUT: SEM.dsrp_adding_genes_trun_stack_bic
        *
        * All of the models pick the baseline as the best model.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_identify_best_model_newlink.sas';

/* SEM Moving genes around Genome-Wide Genes */
    /* Combine BIC from moving all genes */
        * For each gene, add the baseline model BIC. Format dataset
        *
        * libname addgen '!MCLAB/cegs_sem_sd_paper/adding_genes/dspr_gene_level_newlink/sas_data'
        * 
        * INPUT: addgen.*
        *
        * DATASET: SEM.dsrp_gene_level_add_newlink
        *          SEM.dsrp_gene_level_add_newlink_sbs
        * 
        * FILES: !MCLAB/cegs_sem_sd_paper/adding_genes/dsrp_gene_level_add_newlink.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_newlink.sas';

    /* Identify best model for each gene and see what stands out */
        * 
        * INPUT: SEM.dsrp_gene_level_add_newlink
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_identify_best_model_newlink.sas';

