/*******************************************************************************
* Filename: Makefile_splicing.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Look for evidence of splicing effects. Motivation: we are trying
* to understand why InR is loading onto Sxl. One possible explanation is that
* InR is differentially spliced, and Sxl's association with the splicing
* machinery could be why it is connected.
*
*******************************************************************************/

/* Libraries */
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname cegs '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Look at overall splicing in Sex det genes and InR */
    /* Build Table for Analysis */
        * Create a stacked dataset with sex det and InR genes.
        *
        * INPUT:CEGS.ccfus_norm_centered
        *       SEM.cegs_flag_sex_det
        *       DMEL551.symbol2coord
        *       DMEL551.fb551_si_fusions_unique_flagged
        *
        * DATASET: SEM.cegsV_splice_data
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_build_table_to_test_splicing.sas';

    /* Run a basic splice model */
        * Run a basic splice model to determine if there is a splicing effect
        * for each gene.
        *
        * INPUT: SEM.cegsV_splice_data
        * 
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/splicing/*.html
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_differential_splice_analysis.sas';

/* CegsV InR Splicing */
    * Trying to figure out if InR is differentially spliced by genotype.;

    /* Run statisical models */
        * Run a set of linear models to decide which fusions and genotypes to
        * look at.
        *
        * INPUT: SEM.cegsV_splice_data
        *
        * DATASET:
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_InR_splicing_model.sas';


    /* Buld Table showing of InR fusions with normalized value */
        * Build a simple sbs table with InR fusions and their normalized
        * values. 
        *
        * Note: that several InR fusions are missing because they
        * were removed during the normalization process.
        * 
        * INPUT: DMEL551.symbol2coord
        *        DMEL551.fb551_si_fusions_unique_flagged
        *        SEM.cegs_virgin_norm_cent
        *
        * DATASET: SEM.inr_coverage
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_create_cegsV_inr_table.sas';

    /* InR fusions table with fusion level apn counts */
        * Build table of the raw counts for the inr fusions. Add on the overall
        * gene value that I am using for SEMs.
        * 
        * INPUT: DMEL551.fb551_si_fusions
        *        CEGS.flag_drop_fusion_by_line
        *        CEGS.line_norm_centered
        *        SEM.cegsV_by_gene_sbs
        *
        * DATASET: SEM.inr_gene_and_fusion
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/reports/inr_fusions.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_cegsV_inr_raw_table.sas';

    /* Look at Correlation between fusion level apn and gene level */
        * Look at if there is a single fusion that is driving the relationship
        * between InR and Sxl. See how the dsx relationships change.
        *
        * INPUT: SEM.inr_gene_and_fusion
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/splicing/correlations.html
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_cegsV_inr_corr_with_fusion.sas';

    /* Create table of fusion ratios */
        * Make table of fusion ratios. I will use this data to make cell
        * plots. I am specifically interested in the ratio of S57607/S57601
        * and S57607/S57607. 
        *
        * INPUT: SEM.inr_gene_and_fusion
        *
        * DATASET: SEM.inr_fusion_ratio
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_calc_inr_ratios.sas';

/* dsxNull InR Splicing  */
    * Sergey had an idea that we could take the DSXNull data and look and see
    * if there is differential splicing between compared to the wildtype. To do
    * this I need to pull the InR fusions that are present in the cegsV
    * dataset. Then run the differential splicing model from above with control
    * vs Null.

    /* Build DSX table */
        * Build a table containing the InR fusions that are present in the cegsV data.
        *
        * INPUT: SEM.InR_gene_and_fusion
        *        DSX.dsx_all 
        *
        * DATASET: SEM.dsxnull_inr_fusions
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_build_dsx_table_to_test_splicing.sas';

    /* Run a basic splice model */
        * Run a basic splice model to determine if there is a splicing effect
        * between treatments
        *
        * INPUT: SEM.dsxnull_inr_fusions
        * 
        * FILE: 
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/splicing_dsxnull_differential_splice_analysis.sas';
