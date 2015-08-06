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
libname fru '!MCLAB/arbeitman/arbeitman_fru/sasdata';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';
libname genelist '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Import Luo dataset */
    * The Luo2011 dataset contains a list of genes that have dsx binding sites.
    * Significant genes suggest biological evidence using DAMID.
    *
    * INPUT: !MCLAB/useful_dmel_data/gene_lists/Luo2011_Table1_dsx_binding_sites.csv
    *
    * DATASET: GENELIST.Luo2011_dsx_binding_site
    ;
    %include '!MCLAB/useful_dmel_data/gene_lists/sas_programs/import_luo2011_dsx_binding_sites.sas';

/* Import Chang2011 */
    * The Chang2011 dataset contains a list of genes that are affected by tra
    * pseudomales. It also contains tra binding site information.
    *
    * INPUT: !MCLAB/useful_dmel_data/gene_lists/chang2011/Chang2011_TableS10_deg.csv
    *        !MCLAB/useful_dmel_data/gene_lists/chang2011/Chang2011_TableS20_tra_bs.csv
    *
    * DATASET: GENELIST.chang2011_tra
    ;
    %include '!MCLAB/useful_dmel_data/gene_lists/sas_programs/import_chang2011_tra.sas';

/* Import Goldman2007 */
    * The Goldman2007 dataset contains a list of genes that are affected by tra
    * pseudomales. It also contains tra binding site information.
    *
    * INPUT: !MCLAB/useful_dmel_data/gene_lists/Goldman2007/Goldman2007_TableS1_MF_deg.csv
    *        !MCLAB/useful_dmel_data/gene_lists/Goldman2007/Goldman2007_TableS2_genes_ds_tra.csv
    *        !MCLAB/useful_dmel_data/gene_lists/Goldman2007/Goldman2007_TableS5_genes_ds_tra_not_dsx.csv
    *
    * DATASET: GENELIST.Goldman2007_tra
    ;
    %include '!MCLAB/useful_dmel_data/gene_lists/sas_programs/import_goldman2007_tra.sas';

/* Import McIntyre2006 */
    * The McIntyre2006 dataset contains a list of genes that sex bias in
    * splicing.
    *
    * INPUT: !MCLAB/useful_dmel_data/gene_lists/McIntyre2006/McIntyre2006_TableS1.csv
    *
    * DATASET: GENELIST.McIntyre2006
    ;
    %include '!MCLAB/useful_dmel_data/gene_lists/sas_programs/import_mcintyre2006_sex_biased_splicing.sas';

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
    *        SEM.cegsV_ag_yp2_added_genes
    *
    * DATASET: SEM.validation_set (334795 obs)
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_merge_all.sas';

/* Validation Freqs */
    * Now that I have a dataset that contains all of the external dataset, run
    * some freqs to see what is going on.
    *
    * INPUT: SEM.validation_set
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_freqs.sas';

/* Overlap between GRN expansion and DsxNull */
    * Create table of the paper that contains the genes that overlap between
    * the 1390 added genes and the genes from dsxNull.
    *
    * INPUT: SEM.validation_set
    *
    * DATASET:
    *
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_grn_expanstion_dsx_table.sas';

/* Enrichment of Sex Biased Splicing */
    * Sex det pathway is a splicing cascade, so we may expect that genes
    * previously shown as sex bias, may be effected.
    *
    * INPUT: SEM.validation_set
    * 
    * DATASET:
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/validation_enrich_sex_bias.sas';

