/********************************************************************************
* Create induced and repressed gene list for the different treatment types
* based on ANOVA results.
*
* Michelle had done this on her side, but I created all of these to test and
* make sure there were no errors. I compared my output list with the ones she
* sent and everyhing matched up.
********************************************************************************/

libname fru '!MCLAB/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Create FruMA Male Gene Lists */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Average Fold difference >= 2 for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruA_male_fusions_ind_rep
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_male_ind_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_male_rep_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruA_male.sas';

/* Create FruMA Male Gene Lists with no fold change criteria */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Directions of difference is the same for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruA_male_fus_ind_rep_nfold
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_male_ind_1fold_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_male_rep_1fold_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruA_male_nofold.sas';

/* Create FruMA Female Gene Lists */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Average Fold difference >= 2 for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruA_female_fusions_ind_rep
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_female_ind_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_female_rep_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruA_female.sas';

/* Create FruMA Female Gene Lists with no fold change criteria */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Directions of difference is the same for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruA_female_fus_ind_rep_nfold
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_female_ind_1fold_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_A_female_rep_1fold_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruA_female_nofold.sas';

/* Create FruMB Male Gene Lists */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Average Fold difference >= 2 for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruB_male_fusions_ind_rep
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_male_ind_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_male_rep_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruB_male.sas';

/* Create FruMB Male Gene Lists with no fold change criteria */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Directions of difference is the same for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruB_male_fus_ind_rep_nfold
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_male_ind_1fold_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_male_rep_1fold_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruB_male_nofold.sas';

/* Create FruMB Female Gene Lists */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Average Fold difference >= 2 for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruB_female_fusions_ind_rep
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_female_ind_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_female_rep_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruB_female.sas';

/* Create FruMB Female Gene Lists with no fold change criteria */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Directions of difference is the same for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruB_female_fus_ind_rep_nfold
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_female_ind_1fold_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_B_female_rep_1fold_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruB_female_nofold.sas';

/* Create FruMC Male Gene Lists */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Average Fold difference >= 2 for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruC_male_fusions_ind_rep
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_male_ind_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_male_rep_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruC_male.sas';

/* Create FruMC Male Gene Lists with no fold change criteria */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Directions of difference is the same for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruC_male_fus_ind_rep_nfold
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_male_ind_1fold_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_male_rep_1fold_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruC_male_nofold.sas';

/* Create FruMC Female Gene Lists */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Average Fold difference >= 2 for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruC_female_fusions_ind_rep
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_female_ind_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_female_rep_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruC_female.sas';

/* Create FruMC Female Gene Lists with no fold change criteria */
    * Criteria:
    *       FDR <0.2 for all 4 fru null comparisons
    *       Directions of difference is the same for each of the fru comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.fruC_female_fus_ind_rep_nfold
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_female_ind_1fold_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/FruM_C_female_rep_1fold_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_fruC_female_nofold.sas';

/* Create Fru Null Male Gene List */
    * Criteria:
    *       FDR >0.2 for all 4 fru null comparisons
    *       Average Fold difference > 1 for each of the fru null comparisons
    *
    * INPUT: FRU.results_plus_gov2
    * 
    * DATASET: FRU.Null_male_fusions_ind_rep
    *
    * OUTPUT: '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/Induced_Fru_m_null_jmf.tab'
    *         '!MCLAB/arbeitman/arbeitman_fru_network/exported_data_from_michelle/Repressed_Fru_m_null_jmf.tab'
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_null_male.sas';

/* Combine gene lists with Results Table */
    * Combine the results table with the fusion level flags if a fusion was
    * considered induced or repressed.
    *
    * INPUT: FRU.results_by_fusion
    *        FRU.fruA_male_fusions_ind_rep
    *        FRU.fruA_female_fusions_ind_rep
    *        FRU.fruB_male_fusions_ind_rep
    *        FRU.fruB_female_fusions_ind_rep
    *        FRU.fruC_male_fusions_ind_rep
    *        FRU.fruC_female_fusions_ind_rep
    *        FRU.Null_male_fusions_ind_rep
    *
    * DATASET: FRU.results_with_flag_ind_rep
    ;
    %include '!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/DEG_list_combine_results_and_flag_ind_rep.sas';

