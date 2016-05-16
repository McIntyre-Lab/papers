/*******************************************************************************
* Filename: Makefile_100_genome_simulation.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Prepare ASE counts for input in the Bayesian machine using
* emprical q-values.
*******************************************************************************/

/* Libraries */
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Import sam-compare and prep for bayesian machine */
    * Import sam-compare results for all 100 lines
    * 
    * INPUT: !MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/ase_counts_fb551/ase_counts_Line&ii..csv
    *        WORK.design_file
    *
    * DATASET: WORK.all_ase
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/100_genome_simulation_import_sam-compare.sas';

/* Merge Proporation Line */
    * Take the proporation of Line / total from each line and merge into a
    * sbs.
    *
    * INPUT: WORK.all_ase
    *
    * DATASET: WORK.wide
    *          CEGS.fb551_100_genome_flag_line_bias
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/fb551_100_genome_flag_line_bias.csv
    *       !MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/fb551_100_genome_bias_wide.csv
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/100_genome_simulation_merge_proportions.sas';
    
/* Create SAS version of drop list */
    * I made a drop list in python, I now want to import the flags into SAS. See
    * $PROJ/scripts/100_genome_simulation/genome_ambiguity_summary.ipynb for
    * more information.
    *
    * FILE:
    * '!MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/flag_exonic_region_w_and_wo_bias_100_genome_simulation.xls'
    *
    * DATASET: CEGS.exon_drop_list_100_genome
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/100_genome_simulation_create_drop_list.sas';
