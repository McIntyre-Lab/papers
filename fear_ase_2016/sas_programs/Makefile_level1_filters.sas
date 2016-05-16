/*******************************************************************************
* Filename: Makefile_level1_filters.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Run level 1 filters
*
*******************************************************************************/

/* Libraries */
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);

/* Create CEGS genotype list */
    * Create a dataset containing a list of the 68 genotypes for the CEGS ASE
    * paper
    * INPUT: !MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/CEGS_list_68_lines.txt
    * 
    * DATASET: CEGS.genotype_list
    *          CEGS.ase_design_file
    ;
    *%include '!MCLAB/cegs_ase_paper/sas_programs/create_genotype_list.sas';

data design_file;
    set CEGS.genotype_list;
    run;

/* Import SNP VCFs */
    * Import split SNP and Indel vcf files into SAS
    *
    * INPUT: WORK.design_file
    *        !HOME/sandbox/cegs_ase_paper/snps_by_sample/*.vcf
    *
    * DATASET: WORK.*_vcf
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/import_vcf_files_v2.sas';

/* Count Number of snps between w1118 and each line */
    * Merge w1118 with each line and count the number of snps they share. I
    * think this step can be skipped, but need to check.
    *
    * INPUT: WORK.design_file
    *        WORK.*_vcf
    * 
    * DATASET: THUMP.snp_diffs_w1118_2_line_completes 
    ;
    *%include '!MCLAB/cegs_ase_paper/sas_programs/cnt_snp_diffs_between_w1118_and_line_jmf4.sas';

/* Merge SNP VCF file pairwise with w1118 */
    * Merges the VCF files for each w1118-line pair and flags hets, as well as
    * positions (SNPS, but could be the reference base) in both files or in one
    * file and not the other.
    *
    * INPUT: WORK.design_file
    *        WORK.*_vcf
    * 
    * DATASET: WORK.w1118_2_&ID._vcf
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/merge_vcf_files_4_each_w1118_line_pair_jmf2.sas';

/* Create Filtered SNPs for Masking */
    * Uses 7 Cases to flag SNPs between lines and references. Sets filters for
    * masking. Output file for all line to w1118 comparisons. Do not delete
    * THUMP.w1118_2_&ID._vcf files.
    *
    * INPUT: WORK.design_file
    *        WORK.w1118_2_&ID._vcf 
    *
    * DATASET: THUMP.w1118_2_&ID._vcf_flags2
    *
    * FILE: !HOME/sandbox/ase_lvl1_filtered_vcf_files/w11182&ID._MSKD.vcf
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/create_filtered_vcf_w_w11182line_SNPs_jmf2.sas';

/* Create Filtered for counting */
    * Uses 7 Cases to flag SNPs between lines and references. Sets filters for
    * masking. Output file for all line to w1118 comparisons. 
    *
    * INPUT: THUMP.w1118_2_&ID._vcf_flags2
    *
    * FILE: !HOME/sandbox/ase_lvl1_filtered_vcf_files/w1118_w11182&ID._lvl1.vcf
    *       !HOME/sandbox/ase_lvl1_filtered_vcf_files/&ID._w11182&ID._lvl1.vcf
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/create_filtered_vcf_w_w11182line_SNPs_for_counting_jmf2.sas';
