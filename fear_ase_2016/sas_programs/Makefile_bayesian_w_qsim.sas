/********************************************************************************
* SAS protion of the Bayesian Machine
********************************************************************************/

libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname cegsqc '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/qc/sas_data';
libname design '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Import sam-compare and prep for Bayesian machine */
    * Imports sam-compare results (ignoring contaminated lines)
    * 
    * INPUT: !MCLAB/cegs_sergey/pipeline_output/ase_counts_fb551_updated_fusions/ase_counts_&line._w1118&line._&mating_status._&rep..csv
    *        DESIGN.complete_design_by_rep
    *        CEGSQC.flag_contaminated 
    *
    * DATASET: CEGS.ase_counts_for_bayesian
    *
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_bayesian_machine_import_sam-compare.sas';

/* Create flag analyze */
    * A fusion is analyzable if it has at least one ASE read in a rep.
    * 
    * flag_analyze = 1 if there was ASE information for at least one rep
    * 
    * INPUT: CEGS.ase_counts_for_bayesian
    *
    * DATASET: CEGS.bayesian_flag_analyze
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_machine_flag_analyze.sas';

/* WITHOUT BOTH */
    /* Remove Lines without <3 reps and create SBS */
        * Removes line*mv that have less than 3 replicates. Then create sbs dataset
        * with the num_reps.
        * 
        * INPUT: CEGS.ase_counts_for_bayesian
        *
        * DATASET: CEGS.ase_counts_for_bayesian_sbs
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_machine_remove_and_make_sbs.sas';

    /* Build and export datasets for bayesian machine */
        * (1) Packages dataset with flags and counts
        * (2) Creates a list of lines
        * (3) Creates a list of fusions
        * (4) Takes output dataset and makes mated virgin side-by-side for bayeisan machine
        * 
        * INPUT: CEGS.ase_counts_for_bayesian_sbs
        *        CEGS.ase_qsim_line
        *        CEGS.bayesian_flag_analyze
        *
        * DATASET: CEGS.data_for_bayes_sbs_mv
        *          CEGS.bayesian_fusion_list
        *
        * FILE: !MCLAB/cegs_sergey/pipeline_output/ase_bayesian/ase_dataset_sbs_w_qsim_20140628.csv
        *       !MCLAB/cegs_sergey/pipeline_output/ase_bayesian/bayesian_design_fusion_list.csv
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_machine_build_and_export_w_qsim.sas';

    /* Make a test dataset for Luis */
        * Make a test dataset for Luis
        * 
        * INPUT: CEGS.data_for_bayes_sbs_mv
        * 
        * FILE: !MCLAB/cegs_sergey/pipeline_output/ase_bayesian/r101_mv_sbs_w_qsim_test_dataset.csv
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_luis_test_set.sas';

    /* Export by fusion */
        * On the HPC seems to be most efficient to run the bayesian machine by
        * fusion.
        *
        * INPUT: CEGS.data_for_bayes_sbs_mv
        *
        * FILES: /home/jfear/tmp/fusions/*
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_machine_export_fusions_MV_sbs_w_qsim.sas';

/* WITH BOTH */
    /* Remove Lines without <3 reps and create SBS WITH BOTH */
        * Removes line*mv that have less than 3 replicates. Then create sbs dataset
        * with the num_reps.
        * 
        * INPUT: CEGS.ase_counts_for_bayesian
        *
        * DATASET: CEGS.ase_counts_for_bayesian_both_sbs
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_machine_remove_and_make_sbs_v2.sas';

    /* Build and export dataset for bayesian machine BOTH */
        * (1) Packages dataset with flags and counts
        * (2) Takes output dataset and makes mated virgin side-by-side for bayeisan machine
        * 
        * INPUT: CEGS.ase_counts_for_bayesian_both_sbs
        *        CEGS.ase_qsim_both
        *        CEGS.bayesian_flag_analyze
        *
        * DATASET: CEGS.data_for_bayes_both_sbs_mv
        *
        * FILE: !MCLAB/cegs_sergey/pipeline_output/ase_bayesian/ase_dataset_sbs_w_qsim_with_both_20140702.csv
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_machine_build_and_export_w_qsim_v2.sas';

    /* Make a test dataset for Luis with BOTH */
        * Make a test dataset for Luis
        * 
        * INPUT: CEGS.data_for_bayes_both_sbs_mv
        * 
        * FILE: !MCLAB/cegs_sergey/pipeline_output/ase_bayesian/r101_mv_sbs_w_qsim_test_dataset.csv
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_luis_test_set.sas';

    /* Export by fusion with BOTH */
        * On the HPC seems to be most efficient to run the bayesian machine by
        * fusion.
        *
        * INPUT: CEGS.data_for_bayes_both_sbs_mv
        *
        * FILES: /home/jfear/tmp/fusions/*
        ;
        %include '!MCLAB/cegs_sergey/sas_programs/ase_bayesian_machine_export_fusions_MV_sbs_w_qsim_both.sas';



/********************* RUN BAYESIAN MACHINE ON HPC *********************/

/* Import Bayesian machine results */
    * Import the Bayesian machine results for q = 0.4, 0.5, 0.6
    *
    * FILES: !MCLAB/cegs_sergey/pipeline_output/bayesian/bayesian_results_summary.csv
    *
    * DATASET: CEGS.flag_empirical_bayes
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/bayesian_import_summary.sas';

/* Import Bayesian machine results -- extended */
    * Import the Bayesian machine results for q = 0.4, 0.5, 0.6
    *
    * FILES: !MCLAB/cegs_sergey/pipeline_output/bayesian/bayesian_extended_summary.csv
    *
    * DATASET: CEGS.flag_empirical_bayes_extended
    *
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/bayesian_import_summary.sas';

/* Merge on Gene information and look at distriubtion of theta */
    * In order to make sense out of things, lets merge on gene ids and do some
    * distribution plots of mean theta.
    *
    * INPUT: CEGS.flag_empirical_bayes
    *
    * DATASET: CEGS.flag_empirical_bayes_w_symbol
    *
    * RSCRIPT: $MCLAB/cegs_sergey/r_programs/ase_bayesian_plot_mini_empirical_thetas.R
    * 
    * FILE: !MCLAB/cegs_sergey/reports/ase/ase_bayesian_preliminary_results.csv
    *       $MCLAB/cegs_sergey/reports/ase/ase_bayesian_distribution_theta.png
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/bayesian_merge_gene_symbol.sas';

/* Merge APN to Empirical Bayesian Results */
    * Merge on apn so that we can look at the distibution for fusions*line*mv
    * that have AI no matter what the q is.
    *
    * INPUT: CEGS.flag_empirical_bayes
    *        CEGS.ccfus_stack
    *
    * FILES: !MCLAB/cegs_sergey/bayesian_analysis/empirical_bayesian_AI_vs_apn.csv
    ;
    %include '!MCLAB/cegs_sergey/sas_programs/bayesian_merge_apn.sas';
