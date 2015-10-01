/*******************************************************************************
* Filename: Makefile_level2_filters.sas
*
* Author: Alison Gerken | agerken@ufl.edu
*
* Description: Run level 2 filters
*
*******************************************************************************/

/* Libraries */
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname thump '!HOME/thumper/cegs_ase_paper/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

data design_file;
    set CEGS.genotype_list;
    run;

/* Import RNA support */
    * Import the summarized RNA counts for each genotype
    *
    * INPUT: WORK.design_file
    *        !HOME/sandbox/cegs_ase_paper/ase_masked_aln_rna_support/&ID..csv
    *
    * DATASET: WORK.&ID._counts
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/level2_import_rna_support.sas';

/* Merge to Level 1 Filters */
    * Merge RNA counts to level 1 filters and make level 2 filter
    * Lvl2 Filter = RNA CVG > 0 or DNA CVG > 5
    *
    * I am also removing any overlapping variants, because they are ambiguous
    *
    * INPUT: WORK.design_file
    *        THUMP.w1118_2_&ID._vcf_flags2
    *        THUMP.w1118_2_&ID._indel_vcf_flags2
    *        WORK.&ID._counts
    * 
    * DATASET: THUMP.flag_lvl2_w1118_2_&ID.
    *          THUMP.perm_mask_w11182&ID.
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/level2_merge_level1_filters.sas';

/* Split and Export UPD VCFs */
    * Split VCFs for w1118 and line, keeping the SNPs and indels together. Then
    * export a vcf for updating.
    *
    * INPUT: WORK.design_file
    *        THUMP.flag_lvl2_w1118_2_&ID.
    *        THUMP.perm_mask_w11182&ID.
    * 
    * FILE: !HOME/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files/&ID._w11182&ID._UPD.vcf
    *       !HOME/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files/w1118_w11182&ID._UPD.vcf
    *       !HOME/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files/pos_to_permMask_w11182&ID..bed
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/level2_export_UPD_vcf.sas';
