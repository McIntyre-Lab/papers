/*******************************************************************************
* Filename: Makefile_bayesian_emp_dataprep.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Prepare ASE counts for input in the Bayesian machine using
* emprical q-values.
*
*******************************************************************************/

/* Libraries */
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

data design_file;
    set CEGS.ase_design_file;
    run;

/* Import sam-compare and prep for bayesian machine */
    * Import sam-compare results
    * 
    * INPUT: !MCLAB/cegs_ase_paper/pipeline_output/ase_counts_fb551_updated_fusions/ase_counts_&line._&mating_status.&rep..csv
    *        WORK.design_file
    *
    * DATASET: WORK.all_ase
    *
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/emp_bayesian_import_sam-compare.sas';

/* Prepare dataset for bayesian machine */
    * This script does the following:
    *   (1) Remove line*mv that have less than 3 
    *   (2) Create a flag_analyze for each fusion if APN > 0 for at least 1 rep
    *   (3) Package dataset with flags and counts
    *   (4) Export full dataset
    *   (5) Create a list of lines
    *   (6) Create a list of fusions
    * 
    * INPUT: WORK.all_ase
    *
    * DATASET: CEGS.emp_bayesian_input
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/emp_bayesian/input/ase_dataset_for_bayesian.csv
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/emp_bayesian_prep_data_for_machine.sas';

/* Import Bayesian machine results */
    * Import the Bayesian machine results for q = 0.4, 0.5, 0.6
    *
    * FILES: !MCLAB/cegs_ase_paper/pipeline_output/emp_bayesian/PG_model/PG_emp_bayesian_results.csv
    *
    * DATASET: CEGS.emp_bayesian_results
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/emp_bayesian_import_results.sas';

/* Look at distribution of AI calls (Thetas) */
    * Look at the distribution of the thetas to determine how the empirical
    * Bayesian estimates performed. Going to merge on APN and flags from the
    * 100 genome simulation to remove fusions that will not perform well in the
    * bayesian machine.
    *
    * NOTE: right now I have it set to if a fusion showed bias in any number of
    * simulated lines. I may want to tweak this to say 50% of lines;
    *
    * INPUT: CEGS.emp_bayesian_results
    *        CEGS.emp_bayesian_input
    *        CEGS.fb551_100_genome_bias_counts
    *
    * DATASET: CEGS.emp_bayesian_results_w_flags
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/emp_bayesian/PG_model/emp_for_plotting.csv
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/emp_bayesian_merge_gene_symbol_and_apn.sas';
