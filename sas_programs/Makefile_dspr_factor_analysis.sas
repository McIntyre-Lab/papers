libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Factor Analysis */
    /* Look at Communalities */
        * I need to select how many factors to keep. In order to do this I am
        * looking at the communalities to determine where all factors have say
        * >.75. To this I am running many factor analysis and output the
        * communalities. Then I can look at this table and see which genes
        * improve together with different communalities.
        *
        * Somewhere ~14-17 factors looks good.
        *
        * INPUT: SEM.dsrp_sex_det_sbs_combine_sym
        *
        * DATASET: SEM.dsrp_sex_det_factor_communal
        * 
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_communal/dsrp_sex_det_factor_analysis_&num..log
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_communal/dsrp_sex_det_factor_analysis_&num..lst
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_factor_analysis_communalities.sas';

    /* Factor Analysis on All Isoforms*/
        * Perform a simple factor analysis using the isoforms. Looking for
        * factor loadings as well as making sure all of the matrices perform
        * well with the inversions and rotations. Keep 16 factors given the
        * above communalities.
        * 
        * INPUT: SEM.dsrp_sex_det_sbs_combine_sym
        *
        * SCRIPT: $MCLAB/cegs_sem_sd_paper/scripts/plot_fa_grid.py 
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.csv
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.png
        *
        * NOTE: There are some within gene isoforms that are loading on the
        * same factor. It may be worth combining these for simplicity.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_factor_analysis_all_iso.sas';

    /* Factor Analysis on Collapsed Isoforms*/
        * Perform a simple factor analysis using the collapsed isoforms.
        * Looking for factor loadings as well as making sure all of the
        * matrices perform well with the inversions and rotations. Keep 16
        * factors given the above communalities.
        * 
        * Some isoforms always load on the same factor. So I am iterating back
        * and forth with the data preparation to combine these genes. Namely
        * "Spf45" and "sqd" seem to be the biggest problems.
        *
        * INPUT: SEM.dsrp_sex_det_sbs_combine_sym
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.csv
        *       $MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.png
        *
        * NOTE: There are some within gene isoforms that are loading on the
        * same factor. It may be worth combining these for simplicity.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_factor_analysis_combined_iso.sas';

    /* Factor Analysis on Genes*/
        * Perform a simple factor analysis using genes.  Looking for factor
        * loadings as well as making sure all of the matrices perform well with
        * the inversions and rotations. Keep 8 factors. 
        * 
        * INPUT: SEM.dsrp_sex_det_sbs_gene_level_sym
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.csv
        *       $MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.png
        *
        * NOTE: There are some within gene isoforms that are loading on the
        * same factor. It may be worth combining these for simplicity.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/dsrp_sex_det_factor_analysis_genes.sas';
