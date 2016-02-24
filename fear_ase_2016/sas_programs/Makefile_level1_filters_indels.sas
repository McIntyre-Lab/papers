/*******************************************************************************
* Filename: Makefile_level1_filters.sas
*
* Author: Justin M Fear | jfear@ufl.edu ; Alison Gerken | agerken@ufl.edu
*
* Description: Run level 1 filters on INDELS
*
*******************************************************************************/

/* Libraries */
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper/sas_data';
libname DMEL '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Create CEGS genotype list */
    * Create a dataset containing a list of the 68 genotypes for the CEGS ASE
    * paper. Same as SNP step.
    *
    * INPUT: !MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/CEGS_list_68_lines.txt
    * 
    * DATASET: CEGS.genotype_list
    ;
    *%include '!MCLAB/cegs_ase_paper/sas_programs/create_genotype_list.sas';

data design_file;
    set CEGS.genotype_list;
    run;

/* Import VCF */
    * Import split vcf files into SAS. Split previously using splitVcfBySample.py
    *
    * INPUT: CEGS.genotype_list
    *        !HOME/sandbox/cegs_ase_paper/indels_by_sample/*.vcf
    *
    * DATASET: WORK.*_indel_vcf
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/import_indel_vcf_files.sas';

/* Count Number of snps between w1118 and each line */
    * Merge w1118 with each line and count the number of snps they share. I
    * think this step can be skipped, but need to check.
    *
    * INPUT: CEGS.genotype_list
    *        WORK.*_indel_vcf
    * 
    * DATASET: THUMP.indel_diffs_w1118_2_line_completes 
    ;
    *%include '!MCLAB/cegs_ase_paper/sas_programs/cnt_indel_diffs_between_w1118_and_line.sas';

/* Merge VCF file pairwise with w1118 */
    * Merges the VCF files for each w1118-line pair and flags hets, as well as
    * positions (INDELS, but could be the reference base) in both files or in one
    * file and not the other.
    *
    * INPUT: CEGS.genotype_list
    *        WORK.*_indel_vcf
    * 
    * DATASET: WORK.w1118_2_&ID._indel_vcf
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/merge_indel_vcf_files_4_each_w1118_line_pair_jmf2.sas';

/* Create Filtered Indels */
    * Uses 7 Cases to flag SNPs between lines and references. Sets filters for
    * masking. Output file for all line to w1118 comparisons. Do not delete
    * THUMP.w1118_2_&ID._vcf files.
    *
    * INPUT: CEGS.genotype_list
    *        WORK.w1118_2_&ID._indel_vcf
    *
    * DATASET: THUMP.w1118_2_&ID._indel_vcf_flags2
    *
    * FILE: !HOME/sandbox/cegs_ase_paper/ase_lvl1_filtered_vcf_files/w1118_w11182&ID._lvl1_indel.vcf
    *       !HOME/sandbox/cegs_ase_paper/ase_lvl1_filtered_vcf_files/&ID._w11182&ID._lvl1_indel.vcf
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/create_filtered_vcf_w_w11182line_indels_jmf3.sas';
