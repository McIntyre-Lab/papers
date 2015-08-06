/*******************************************************************************
* Filename: Makefile_cegs_validation.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';
libname genelist '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Merge Valdiation datset to Adding Genes */
    * There are several lines of external data that we are using for
    * validation. Merge all of these datasets onto the 1390 genes added by the
    * SEM.
    *
    * INPUT: GENELIST.Luo2011_dsx_binding_site
    *        GENELIST.chang2011_tra
    *        GENELIST.Goldman2007_tra
    *        GENELIST.McIntyre2006
    *        DMEL551.symbol2fbgn
    *        DMEL530.symbol2fbgn
    *        SEM.dsxNullf_induced
    *        SEM.dsxNullf_repressed
    *        SEM.ctrl_female_induced
    *        SEM.ctrl_female_repressed
    *        SEM.cegsV_ag_yp2_stack_bic
    *        SEM.cegsV_ag_model_design_file
    *        SEM.cegsV_gene2fbgn
    *        SEM.cegsV_ag_w_bic_cutoff
    *
    * DATASET: SEM.validation_set_BIC12 (334780 obs)
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_merge_all_BIC12.sas';

/* Validation Freqs */
    * Now that I have a dataset that contains all of the external dataset, run
    * some freqs to see what is going on.
    *
    * INPUT: SEM.validation_set_BIC12
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_freqs_BIC12.sas';

/* DSX Validation Freqs */
    * Focusing on DSX models, what are the FREQS for BS and all things dsx.
    *
    * INPUT: SEM.validation_set_BIC12
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_freqs_BIC12.sas';

/* Overlap between GRN expansion and DsxNull */
    * Create table of the paper that contains the genes that overlap between
    * the 754 added genes and the genes from dsxNull.
    *
    * INPUT: SEM.validation_set_BIC12
    *
    * DATASET:
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_grn_expanstion_dsx_table_BIC12.sas';

/* Enrichment of Sex Biased Splicing */
    * Sex det pathway is a splicing cascade, so we may expect that genes
    * previously shown as sex bias, may be effected.
    *
    * INPUT: SEM.validation_set_BIC12
    * 
    * DATASET:
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_enrich_sex_bias_BIC12.sas';
