libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Full Covariance Model with Yp2 */
    * All exogenous variables are allowed to covary
    *
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_estimates.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_model.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/full_cov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_sem_full_cov.sas';

/* Full Covariance Model with Yp1 instead of Yp2 */
    * All exogenous variables are allowed to covary
    *
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_estimates_yp1.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_model_yp1.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/full_cov_model_yp1.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_sem_full_cov_yp1.sas';

/* Full Covariance Model with Yp3 instead of Yp2 */
    * All exogenous variables are allowed to covary
    *
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_estimates_yp3.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/full_cov_model_yp3.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/full_cov_model_yp3.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_sem_full_cov_yp3.sas';


/* Gene Covariance Model with Yp2 */
    * covariances between genes constrained to 0
    *
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/gene_cov_estimates.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/gene_cov_model.lst
    *                  !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/cegsV_gene_cov_raw_residuals.csv
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/gene_cov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_sem_gene_cov.sas';

/* Modified No Covariance Model with Yp2 */
    * covariances between genes constrained to 0
    *
    * The modified Gene covariance matrix is allowing the splicing factors snf,
    * Spf45, fl(2)d to covary with tra2.
    *
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/modified_gene_cov_estimates.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/modified_gene_cov_model.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/modified_gene_cov_model.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_sem_modified_gene_cov.sas';

/* No Covariance Model with Yp1 instead of Yp2 */
    * covariances between genes constrained to 0
    *
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/gene_cov_estimates_yp1.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/gene_cov_model_yp1.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/gene_cov_model_yp1.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_sem_gene_cov_yp1.sas';

/* No Covariance Model with Yp3 instead of Yp2 */
    * covariances between genes constrained to 0
    *
    * INPUT: SEM.cegsV_by_gene_sbs
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/gene_cov_estimates_yp3.csv
    *
    * OUTPUT SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/gene_cov_model_yp3.lst
    *
    * LOGS SAVED AS: !MCLAB/cegs_sem_sd_paper/analysis_output/sem/cegsV/logs/gene_cov_model_yp3.log
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_sem_gene_cov_yp3.sas';
