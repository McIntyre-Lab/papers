/* Use RANK */
    * Want to test and see if expression level is driving the partial
    * correlation results. To look into this I will use rank instead of
    * expression level.
    ;

    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);

    /* Create Rank Datasets */
        * Calculate RANKS for each sample using the entire dataset and the Sex
        * Det Subset. 
        *
        * INPUT: SEM.dsrp_stack
        *
        * DATASET: SEM.dsrp_stack_rank
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsrp_create_ranks.sas';

    /* Re-create Sex det from Ranks */
        * Using the Full rank dataset I want to create a sex det dataset to see
        * if it matters. Instead of breaking the FA and combining genes, I am
        * just going to do it all there.
        *
        * INPUT: SEM.dsrp_stack_rank
        *
        * DATASET: SEM.dsrp_sexdet_stack_from_rank
        *          SEM.dsrp_sexdet_sbs_sym_from_rank
        *          SEM.dsrp_sexdet_comb_sym_from_rank
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsrp_create_sex_det_from_ranks.sas';

    /* Rank SEM */
        * SEM analysis of the Sex Determination pathway. There are three
        * different SEMs, (1) un-constrained: covariances between exogenous
        * variables estimated, (2) constrained: covarainces between genes
        * constrained to 0, (3) partially constrained: non-significant (tCrit
        * >=1.96) covarainces from [1] constrained to 0.
        *
        * INPUT: SEM.dsrp_sexdet_sbs_comb_sym_from_rank
        *
        * FILES: !MCLAB/cegs_sem_sd_paper/rank_sem_output/unconstrained_estimates.csv
        *        !MCLAB/cegs_sem_sd_paper/rank_sem_output/constrained_estimates.csv
        *        !MCLAB/cegs_sem_sd_paper/rank_sem_output/partially_constrained_estimates.csv
        *        !MCLAB/cegs_sem_sd_paper/rank_sem_output/spf_isoforms_constrained_estimates.csv
        *        !MCLAB/cegs_sem_sd_paper/rank_sem_output/constrained_w_spf_estimates.csv
        *
        * OUTPUT SAVED AS: !MCLAB/rank_sem_gof/rank_sem_output/unconstrained_model.lst
        *                  !MCLAB/rank_sem_gof/rank_sem_output/constrained_model.lst
        *                  !MCLAB/rank_sem_gof/rank_sem_output/partially-constrained_model.lst
        *                  !MCLAB/rank_sem_gof/rank_sem_output/spf_isoforms_constrained_model.lst
        *                  !MCLAB/rank_sem_gof/rank_sem_output/constrained_w_spf_model.lst
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsrp_sex_det_rank_sem.sas';

    /* RANK GGM on all genes is DSRP */
        * Create and export dataset for gaussian graphical models using the sex
        * determination subset. Then run GGM R script. GGM allows you to select
        * edges by an FDR cutoff or by specifiying a number of edges. Two
        * different FDR cutoffs were used (0.2, 0.8) along with outputing the
        * top 20 edges. GGM was run on the entire DSRP dataset.
        * 
        * INPUT: SEM.rank_dsrp_stack
        *
        * FILES: $MCLAB/cegs_sem_sd_paper/analysis_output/rank_ggm/dsrp_ggm_isoforms_FDR2.dot
        *
        * NOTE:  I could not get sas to export all (11065) columns. I had to
        * export the intermediate csv file using JMP genomics.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsrp_create_rank_ggm_dataset_DSRP.sas';

    /* Gaussian Graphic Model on sex det subset */
        * Create and export dataset for gaussian graphical models using the sex
        * determination subset. Then run GGM R script. GGM allows you to select
        * edges by an FDR cutoff or by specifiying a number of edges. Two
        * different FDR cutoffs were used (0.2, 0.8) along with outputing the
        * top 20 edges. Also Isoforms and gene level variables are used.
        *
        * INPUT: SEM.dsrp_sexdet_comb_sym_from_rank
        * 
        * RSCRIPTS: !MCLAB/cegs_sem_sd_paper/r_programs/dsrp_rank_ggm_sex_det_subset.R
        *
        * FILES: $MCLAB/cegs_sem_sd_paper/analysis_output/rank_ggm/dsrp_ggm_sex_det_isoforms_FDR2.dot
        *        $MCLAB/cegs_sem_sd_paper/analysis_output/rank_ggm/dsrp_ggm_sex_det_isoforms_FDR2.png 
        *        $MCLAB/cegs_sem_sd_paper/analysis_output/rank_ggm/dsrp_ggm_sex_det_isoforms_TOP20.dot
        *        $MCLAB/cegs_sem_sd_paper/analysis_output/rank_ggm/dsrp_ggm_sex_det_isoforms_TOP20.png 
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsrp_create_rank_ggm_dataset_sex_det_subset.sas';

/* DSX Null Dataset */
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel548 '!MCLAB/useful_dmel_data/flybase548/sas_data';
    libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);


    /* Create List of CEGS Adding Genes */
        * The Adding genes process takes a single gene and then tries it in all
        * possible locations in the network. Create a list of genes that had at
        * least one model that was better than the baseline.
        *
        * INPUT: SEM.cegsV_ag_nofilter_yp2_sbs;

        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsx_factor_analysis_combined_iso.sas';

    /* Make DSX gene list dataset (dsxNullF repressed) */
        * Using the gene list from the DSXnull experiment, create a side-by-side
        * dataset with those genes.
        *
        * INPUT: SEM.dsxNullf_repressed
        *        DMEL530.symbol2cg
        *        DMEL548.symbol2proteincg
        * 
        * DATASET: SEM.dsrp_dsx_sbs_symbol
        *          SEM.dsxNullf_repressed_cgnumber
        *
        * NOTES: There were a few changes in annotation between 5.30 and 5.48 that had to
        *        be fixed by hand. See code.
        *
        *        After everything was said and done, only 35 genes from the
        *        DSXnullF repressed dataset were found in the dsrp dataset.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsx_create_sbs.sas';

    /* Calculate correlation between DSX gene list and combine */
        * Multicollinearity is going to be a large problem. It will be best to
        * collapse isoforms if they have a large correlation coefficient (>=0.75).
        * 
        * Uses a python script to calculate the correlation between isoforms
        * and combine them (AVERAGE) if they have a correlation >=0.75.
        *
        * INPUT: SEM.dsrp_dsx_sbs_symbol
        *        SEM.dsxNullf_repressed_cgnumber
        *
        * SCRIPTS: $MCLAB/cegs_sem_sd_paper/scripts/clusterByCorr.py 
        *
        * DATASET: SEM.dsrp_dsx_sbs_combine_sym
        *          
        * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform._dsx.csv
        *        !MCLAB/cegs_sem_sd_paper/analysis_output/correlation/&isoform._dsx.log
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsx_calculate_correlation_and_combine.sas';

    /* DSX Factor Communalities */
        * Using the DSX datasets, I have created a gene list with genes repressed
        * when dsxF is removed (dsxNull females). I iterated over 1-40 factors
        * and created a table to communalities. Looking at this table, 25-30
        * factors seems to explain >65% of varaince for all genes.
        *
        * INPUT: SEM.dsrp_dsx_sbs_symbol
        *        SEM.dsrp_dsx_sbs_combine_sym
        * 
        * DATASET: SEM.dsx_factor_all_communal
        *          SEM.dsx_factor_combine_communal;
        * 
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_all_&num..log
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_all_&num..lst
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_combined_&num..log
        *       !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsx_factor_analysis_combined_&num..lst
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsx_factor_communalities.sas';

    /* DSX Factor Analysis */
        * Perform a simple factor analysis using the collapsed isoforms.
        * Looking for factor loadings as well as making sure all of the
        * matrices perform well with the inversions and rotations. Keep 31
        * factors given the above communalities.
        * 
        * INPUT: SEM.dsrp_dsx_sbs_combine_sym
        *
        * FILE: !MCLAB/cegs_sem_sd_paper/analysis_output/factor/dsrp_dsx_factor_analysis_combined_iso.csv
        *
        * NOTE: There are some within gene isoforms that are loading on the
        * same factor. It may be worth combining these for simplicity.
        ;
        %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_dsx_factor_analysis_combined_iso.sas';

/* Overall Gene Expression */
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
    * Figure out how InR and the sex determination gene expression compares
    * with all of the other genes.
    *
    * INPUT: SEM.cegs_by_gene_stack
    *
    * DATASET: /home/jfear/tmp/for_plot.sas7bdat
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_gene_expression.sas';

/* Check CegsV Replicate Effects */
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname norm '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
    * There is concern that the patterns we are seeing may be driven by outlier
    * samples due to the fact that we do not always have 3 biological
    * replicates.
    *
    * INPUT: NORM.ccfus_norm_centered
    *        SEM.cegs_virgin_norm_cent
    *
    * DATASET: SEM.cegsV_line_rep_number
    * 
    * FILE: !MCLAB/cegs_sem_sd_paper/reports/mean_exp_by_rep.pdf
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_test_check_cegsV_rep_effect.sas';

/* Simultaneously add splicing factors */
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
    * I want to know if I can reach saturation to the F1-hybrid dataset. One of
    * our arguments why DSPR did not add any genes is that there was not enough
    * variation, therefore the model was saturated. I am going to try adding
    * all of the splicing factors simultaneously to the model.
    *
    * INPUT: SEM.cegsV_by_gene_sbs 
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/splicing_corr_analysis.lst
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/splicing_factor_analysis.lst
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_sxl_mmc.lst
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_sxl_factor.lst
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_dsx_mmc.lst
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/saturation_ds_sxl_or_dsx_mmc.lst
    *        /home/jfear/sandbox/splice.csv
    *        /home/jfear/sandbox/dsx.csv
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_tests_simultaneously_add_splicing_factors.sas';

/* Simultaneously add nonA factors */
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
    * I want to know what happens if I add all of the genes associated with
    * nonA and heph
    *
    * INPUT: SEM.cegsV_by_gene_sbs 
    *
    * FILES: !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/splicing_corr_analysis.lst
    *        !MCLAB/cegs_sem_sd_paper/analysis_output/saturation/splicing_factor_analysis.lst
    ;
    %include '!MCLAB/cegs_sem_sd_paper/sas_programs/misc_tests_simultaneously_add_splicing_factors.sas';
