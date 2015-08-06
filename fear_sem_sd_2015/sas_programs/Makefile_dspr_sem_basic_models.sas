libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Full Covariance Model */
    * Basic SEM analysis of the Sex Determination pathway. 
    *
    * Full Covaraince: covariances between exogenous variables estimated
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_estimates.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_model.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/full_cov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_full_cov.sas';

/* Full Covariance Model with Yp1 */
    * Basic SEM analysis of the Sex Determination pathway. 
    *
    * Full Covaraince: covariances between exogenous variables estimated
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_estimates_yp1.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_model_yp1.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/full_cov_model_yp1.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_full_cov_yp1.sas';

/* Full Covariance Model with Yp3 */
    * Basic SEM analysis of the Sex Determination pathway. 
    *
    * Full Covaraince: covariances between exogenous variables estimated
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_estimates_yp3.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_model_yp3.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/full_cov_model_yp3.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_full_cov_yp3.sas';




/* Gene Covariance Model */
    * Basic SEM analysis of the Sex Determination pathway. 
    * Gene covariance: covarainces between genes constrained to 0
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_residuals.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_model.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_gene_cov.sas';

/* Modified Gene Covariance Model */
    * Basic SEM analysis of the Sex Determination pathway. 
    * Gene covariance: covarainces between genes constrained to 0
    *
    * This model allows snf, Spf45, and fl(2)d to covary with tar2
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/modified_gene_cov_estimates.csv
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/modified_gene_cov_residuals.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/modified_gene_cov_model.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/modified_gene_cov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_modified_gene_level_sem_sex_det_gene_cov.sas';

/* Gene Covariance Model with Yp1 */
    * Basic SEM analysis of the Sex Determination pathway. 
    * Gene covariance: covarainces between genes constrained to 0
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_estimates_yp1.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_model_yp1.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_model_yp1.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_gene_cov_yp1.sas';

/* Gene Covariance Model with Yp3 */
    * Basic SEM analysis of the Sex Determination pathway. 
    * Gene covariance: covarainces between genes constrained to 0
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_estimates_yp3.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_model_yp3.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_model_yp3.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_gene_cov_yp3.sas';



/* Truncated Full Covariance Model */
    * Full Covaraince: covariances between exogenous variables estimated without dsx branch of the pathway
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/truncated_fullcov_estimates.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/truncated_fullcov_model.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/truncated_fullcov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_truncated_fullcov.sas';

/* Truncated Gene Covariance Model */
    * Gene covariance: covarainces between genes constrained to 0 without dsx branch of the pathway
    *
    * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym;
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/truncated_gene_cov_estimates.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/truncated_gene_cov_model.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/truncated_gene_cov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_gene_level_sem_sex_det_truncated_gene_cov.sas';




/********************************* DEPRECATED *********************************/
    /* Basic SEM */
        * Basic SEM analysis of the Sex Determination pathway. There are three
        * different SEMs, (1) Full Covaraince: covariances between exogenous
        * variables estimated, (2) Gene covariance: covarainces between genes
        * constrained to 0, (3) no covariance: all covarainces constrained to 0
        *
        * INPUT: SEM.dsrp_sex_det_sbs_combine_sym;
        *
        * DATASET: SEM.dsrp_sex_det_gene_cov_model_bic
        *
        * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_estimates.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_estimates.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/fl2d_isoforms_constrained_estimates.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_w_fl2d_estimates.csv
        *
        * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/full_cov_model.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/gene_cov_model.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/fl2d_isoforms_constrained_model.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_w_fl2d_model.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/truncated_gene_cov_model.lst
        *
        * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/full_cov_model.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/gene_cov_model.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/fl2d_isoforms_constrained_model.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_w_fl2d_model.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/truncated_gene_cov_model.log
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_full_cov.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_gene_cov.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_no_cov.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_fl2d_no_cov.sas';
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_truncated_gene_cov.sas';

    /* SEM check Sxl */
        * Sxl has been grouped into 3 isoforms. Sxl_B and C are annotated as
        * male only isoforms. Sense these are female head, I wanted to look at
        * see if I should include Sxl_B,C in the model. 
        *
        * Found that in general, when I model only SxlA the fit is slightly
        * better. So I will only worry about SxlA for now.
        *
        * INPUT: SEM.dsrp_sex_det_sbs_combine_sym;
        *
        * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_all_sxl.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlB.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlC.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlAB.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlAC.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_estimates_sxlBC.csv
        *
        * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_all_sxl.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlB.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlC.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlAB.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlAC.lst
        *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/constrained_model_sxlBC.lst
        *
        * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_model_all_sxl.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_model_sxlB.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_model_sxlC.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_model_sxlAB.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_model_sxlAC.log
        *                !MCLAB/cegs_sem_sd_paper/analysis_output/sem/logs/constrained_model_sxlBC.log
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_sem_all_sxl.sas';
