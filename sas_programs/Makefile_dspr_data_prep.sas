libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Import DSRP Data */
    * I am using processed microarray data from the DSRP population from King
    * et al. 2014. This data has already been normalized and I will just be
    * modeling these data. The first step is to import the expression matrix.
    *
    * Because the dataset is in the side-by-side and contains >11000 columns,
    * sas was not importing it correctly. To get around this I imported the
    * file in JMP and saved it as a SAS dataset.
    *
    * INPUT: !MCLAB/cegs_sem_sd_paper/original_data/FemaleHeadExpression.txt
    *
    * DATASET: SEM.dsrp;

/* Data Preparation */
    /* Create Sex Determination Subset */
        * Merge on gene annotations and create a subset of only the sex
        * determination genes.
        *
        * INPUT: SEM.dsrp
        *        DMEL548.symbol2proteincg
        *
        * DATASET: SEM.dsrp_sex_det_sbs_cg
        *          SEM.dsrp_sex_det_sbs_symbol
        *          SEM.dsrp_sex_det_stack
        *          SEM.dsrp_stack
        *          SEM.dsrp_sbs_symbol
        *
        * Genes Present In Sex Det Dataset: B52 fl(2)d fru her ix mub ps Psi
        *                                   Rbp1 Rm62 snf Spf45 sqd Sxl tra 
        *                                   tra2 vir Yp1 Yp2 Yp3
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_create_sex_det_dataset.sas';

    /* Create gene list of entire DSRP dataset */
        * It will be useful to have a gene list;
        *
        * INPUT: SEM.dsrp_sbs_symbol;
        *
        * DATASET: SEM.dsrp_gene_list;
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_generate_gene_list.sas';

    /* Factor Analysis on All Isoforms */
        * Perform a simple factor analysis using all of the genes in the sex det
        * pathway.
        *
        * INPUT: SEM.dsrp_sex_det_sbs_symbol
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_all_iso.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_factor_analysis.sas';

    /* Calculate correlation between isoforms in sex det and combine */
        * Multicollinearity is going to be a large problem. It will be best to
        * collapse isoforms if they have a large correlation coefficient (>=0.75).
        * 
        * Uses a python script to calculate the correlation between isoforms
        * and combine them (AVERAGE) if they have a correlation >=0.75.
        *
        * INPUT: SEM.dsrp_sex_det_sbs_symbol
        *
        * SCRIPTS: $MCLAB/cegs_sem_sd_paper/scripts/clusterByCorr.py 
        *
        * DATASET: SEM.dsrp_sex_det_sbs_combine_sym
        *          
        * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform..csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform..log
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform._iso_table.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_calculate_correlation_and_combine.sas';

    /* Calculate correlation between isoforms in all genes and combine */
        * Multicollinearity is going to be a large problem. It will be best to
        * collapse isoforms if they have a large correlation coefficient (>=0.75).
        * 
        * Uses a python script to calculate the correlation between isoforms
        * and combine them (AVERAGE) if they have a correlation >=0.75.
        *
        * INPUT: SEM.dsrp_sbs_symbol
        *
        * SCRIPTS: $MCLAB/cegs_sem_sd_paper/scripts/clusterByCorr.py 
        *
        * DATASET: SEM.dsrp_sbs_combine_sym
        *          
        * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform..csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform._iso_table.csv
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_calculate_correlation_and_combine_v3.sas';

    /* Combine Genome to Gene Level (Average) */
        * Collapse all isoforms to the gene level, so that I can compare it to
        * the cegs dataset.
        *
        * INPUT: SEM.dsrp_sbs_symbol 
        *
        * DATASET: SEM.dsrp_sbs_gene_level_sym
        *          SEM.dsrp_stack_gene_level_sym
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_combine_all_isoforms.sas';

    /* Combine Sex Det to Gene Level (Average) */
        * I will also collapse everything to the gene level. 
        *
        * INPUT: SEM.dsrp_sex_det_sbs_symbol
        *
        * DATASET: SEM.dsrp_sex_det_sbs_gene_level_sym
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_combine_gene_level.sas';

    /* Compare Yp Expression */
        * Look at correlation of gene expression among Yp genes;
        *
        * INPUT: SEM.dsrp_sex_det_sbs_combine_sym 
        * 
        * FILES: !MCLAB/cegs_sem_sd_paper/reports/yp_gene_expression.png
        *
        * NOTE: Creates a corr matrix in the output which needs to be looked at.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_corr_yp.sas';

    /* Compare Sxl Expression */
        * Look at correlation of gene expression among Sxl groups
        *
        * INPUT: SEM.dsrp_sex_det_sbs_combine_sym 
        * 
        * FILES: !MCLAB/cegs_sem_sd_paper/reports/sxl_gene_expression.png
        *
        * NOTE: Creates a corr matrix in the output which needs to be looked at.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_corr_sxl.sas';

