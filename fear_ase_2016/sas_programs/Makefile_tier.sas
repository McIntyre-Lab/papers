/*******************************************************************************
* Filename: Makefile_tier.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Run Type I Error testing for intraspecific simulations.
*
*******************************************************************************/

/* Libraries */
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Intraspecific Data Prep */
    /* Import the Sam-compare RNA counts from virgin */
        * Going to use RNA reads form virgin for simulating data for doing TIER
        * analysis.
        *
        * INPUT: CEGS.ase_design_file
        *
        * FILE: !MCLAB/cegs_ase_paper/pipeline_output/ase_counts_fb551_updated_fusions/ase_counts_&line._&ms.&rep..csv
        *
        * DATASET: WORK.rna_cnts
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_import_rna_cnts_virgin.sas';

    /* Import the Sam-compare DNA counts from read simulation */
        * Going to use DNA counts form the read simulation for simulating data
        * for doing TIER
        *
        * INPUT: CEGS.ase_design_file
        *
        * FILE: !MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/ase_counts_fb551/ase_counts_&line..csv
        *
        * DATASET: WORK.dna_cnts
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_import_simulated_dna_cnts.sas';

    /* Merge and export cnts for TIER simulation */
        * Combine the RNA and simulated DNA counts for and export them for the
        * TIER simulation.
        *
        * INPUT: WORK.rna_cnts
        *        WORK.dna_cnts
        *
        * FILE: !MCLAB/cegs_ase_paper/pipeline_output/typeI_error/input/&line._RNA_sim_DNA_cnts.csv
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_merge_export.sas';

/* Intraspecific TIER Summary */
    /* Import Intraspecific Results */
        * Using intraspecific data simulated from r101, I ran the Bayesian
        * analysis.
        *
        * FILE: !MCLAB/cegs_ase_paper/pipeline_output/typeI_error/output/r101_simulated_NOAI_NOBIASRequal1Poisson_results_summary.csv
        *
        *
        * DATASET: WORK.noai_nobias
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_import_intraspecific.sas';

    /* Calculate TIER percentages for noAI noBIAS */
        * Like the Luis paper, I want to create a table with the TIER percentages
        * for the NoBias NoAI simulated data.
        *
        * INPUT: WORK.noai_nobias
        *
        * Makes FREQS table
        *
        *   Method for AI     TIER (%)
        *   --------------------------
        *   binomial_test     3.66
        *   NB_random_DNA     0.83
        *   PG_q_random_DNA   0.65
        *   PG_q4             27.78
        *   PG_q5             0.46
        *   PG_q6             25.11
        *   flag_all          0.18
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_r101_null_null_table.sas';

/* Intraspecific Misspecification */
    /* Import Bayesian Results */
        * I simulated datasets with differing levels of bias for r101. I ran the Bayesian
        * machine will different amounts of misspecification. Import the results.
        *
        * INPUT: !MCLAB/cegs_ase_paper/pipeline_output/typeI_error/&line._misspecification_results_sim_summary.csv
        *
        * DATASET: CEGS.r361_misspecification
        *          CEGS.r332_misspecification
        *          CEGS.r365_misspecification
        *          CEGS.r101_misspecification
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_import_misspecification.sas';

    /* Intraspecific Misspecification Summary */
        * This script calculates the TIER and formats the table for plotting
        *
        * INPUT: CEGS.&line._misspecification
        *
        * OUTFILE:
        * !MCLAB/cegs_ase_paper/pipeline_output/typeI_error/&line._misspecification_bias_for_plotting.csv
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_misspecification_bias.sas';

    /* Intraspecific Misspecification Binomial Comparison */
        * This script calculates TIER and formats a table for comparing emp
        * Bayes vs Binomial.
        *
        * INPUT: CEGS.&line._misspecification
        *
        * OUTFILE:
        * !MCLAB/cegs_ase_paper/pipeline_output/typeI_error/&line._misspecification_tier_binomial_vs_pg_plotting.csv
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/tier_misspecification_binomial_vs_pg.sas';
