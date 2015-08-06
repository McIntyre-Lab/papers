/*******************************************************************************
* Filename: Makefile_cegs_enrichment.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Run GO enrichment and other enrichment steps if necessary
*
*******************************************************************************/

/* Libraries */
    libname SEM '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
    libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);

/* Prepare data for GO enrichment */
    * I will do go enrichment in JMP Genomics, but need to prep the dataset
    * here.
    *
    * INPUT: SEM.cegsV_ag_w_flags_bic12
    *        DMEL557. genes2go_nodups
    *
    * DATASET: SEM.cegsV_ag_w_flags_bic12_go
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/enrichment_go_data_prep.sas';

/* Chromosomal Enrichment */
    * Are genes added to the network enriched on the different chromosomes?
    *
    * INPUT: SEM.cegsV_ag_w_flags_bic12
    *
    * Did not see any enrichment on for any chromosomes
    * 
    * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/enrichment/chromosome_enrichment_tests.csv
    *       !MCLAB/cegs_sem_sd_paper/analysis_output/enrichment/chromosome_enrichment_tables.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/enrichment_ag_chromosome_enrichment.sas';
    

/* List Enrichment */
    * Using a variety of gene lists test for enrichment.
    *
    * INPUT: SEM.validation_set_bic12
    *
    * DATASET: 
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/enrichment_list_enrichment.sas';
    


