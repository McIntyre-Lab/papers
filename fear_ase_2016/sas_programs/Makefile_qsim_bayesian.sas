/*******************************************************************************
* Filename: Makefile_qsim.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Reads were simulated from the updated genomes. Reads were then
* combined with w1118 to form simulated hybrids. Reads were then aligned to
* both refrences (line, w1118) and sam-compare was run. This set of scripts
* takes the ASE counts and imports them into sas and creates 
*
* qsim = (line / ASE_total)
*
* NOTE: for use in the bayesian machine, q = 1-qsim, because the machine
* calculates theta as (w1118 / ASE_total).
*
*******************************************************************************/

/* Libraries */
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* QSIM from CEGS */
    * For the CEGS project I need to create qsim for each line*fusion_id
    ;
    data design_file;
        set CEGS.genotype_list;
        run;

    /* Import Line Simulation sam-compare and calcualte qsim */
        * This script does the following:
        *   (1) Imports sam-compare results from simulation from lines
        *   (2) calcualtes qsim_tester = tester_total / ASE_total
        * 
        * INPUT: !MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/ase_counts_fb551/ase_counts_&line..csv
        *        design_file
        *
        * DATASET: CEGS.ase_qsim_tester
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/qsim_import_line_simulation_ase_counts.sas';

    /* Check Bias across lines */
        * Do the same genomic regions show bias across all lines.
        *
        * INPUT: CEGS.ase_qsim_tester
        *
        * FILE: !MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/qsim_bias_wide.csv
        *       !MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/num_lines_bias_freqs.rtf
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/qsim_check_bias_across_lines.sas';

    /* Merge QSIM_tester and output for Bayesian Machine */
        * Using the dataset that was prepared for the empricial Bayes, I want
        * to merge on qsim and use this as q in the PG model.
        *
        * INPUT: CEGS.ase_qsim_tester
        *        CEGS.emb_bayesian_input
        *
        * FILE: !MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/input/ase_dataset_for_bayesian_w_qsim.csv
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/qsim_bayesian_prep_data_for_machine.sas';

    /* Import PG qsim results and merge onto the emprical Bayesian results */
        * Import the PG model results using q = qsim. 
        *
        * NOTE: that Luis' program output NA for missing values, which SAS does
        * not like. But since they are missing I am not worrying about them.
        *
        * INPUT: !MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/output/ase_dataset_for_bayesian_w_qsim_summary.csv
        *
        * DATASET: CEGS.qsim_bayesian_results
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/qsim_import_pg_results.sas';
        
    /* Merge qSIM onto emp bayes results */
        * Combine Bayesian results from qsim and emp.
        *
        * When qsim != 1 and qsim != 0 and abs(qsim - 0.5) > 0.05, then use qsim model
        * to estimate AI (flag_AI_combined = 1)
        *
        * Otherwise require all three empirical models to call AI for (flag_AI_combined = 1).
        *
        * INPUT: CEGS.emp_bayesian_results
        *        CEGS.qsim_bayesian_results
        *        CEGS.emp_bayesian_input
        *        CEGS.ase_qsim_tester
        *
        * DATASET: CEGS.qsim_emp_theta_w_flag
        * FILE: !MCLAB/cegs_ase_paper/pipeline_output/ase_bayesian_qsim_emp_theta_w_flag.csv
        ;
        %include '!MCLAB/cegs_ase_paper/sas_programs/qsim_merge_qsim_emp_bayes.sas';
