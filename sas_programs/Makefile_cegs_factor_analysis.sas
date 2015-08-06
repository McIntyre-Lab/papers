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
        * Somewhere ~17-19 factors looks good.
        *
        * INPUT: SEM.cegsV_by_fusion_sex_det_sbs
        *
        * DATASET: SEM.cegsV_sex_det_factor_commun
        * 
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_communal/cegsV_sex_det_factor_analysis_&num..log
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_communal/cegsV_sex_det_factor_analysis_&num..lst
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_factor_analysis_communalities.sas';

    /* Factor Analysis on All Fusions */
        * Perform a simple factor analysis using the isoforms. Looking for
        * factor loadings as well as making sure all of the matrices perform
        * well with the inversions and rotations. Keep 17 factors given the
        * above communalities.
        * 
        * INPUT: SEM.cegsV_by_fusion_sex_det_sbs
        *
        * SCRIPT: $MCLAB/cegs_sem_sd_paper/scripts/plot_fa_grid.py 
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_sex_det_factor_analysis_fusions.csv
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_sex_det_factor_analysis_fusions.png
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_factor_analysis_fusions.sas';

    /* Factor Analysis on Gene Level */
        * Perform a simple factor analysis using genes.  Looking for factor
        * loadings as well as making sure all of the matrices perform well with
        * the inversions and rotations. Keep 13 factors given the scree plot.
        * 
        * INPUT: SEM.cegsV_by_gene_sex_det_sbs
        *
        * SCRIPT: $MCLAB/cegs_sem_sd_paper/scripts/plot_fa_grid.py 
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_sex_det_factor_analysis_gene.csv
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/cegsV_sex_det_factor_analysis_gene.png
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/cegsV_sex_det_factor_analysis_gene.sas';

