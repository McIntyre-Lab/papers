/*******************************************************************************
* Filename: Makefile_ase_summarize.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Pull of the results form empricial Bayes and qsim and estimated
* the amount of ASE. Merge on information about Sex Hierarchy and see what we
* see there.
*
*******************************************************************************/

/* Libraries */
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Filter ASE results */
    * I have explored the ASE results and have found some concerns. After a lot
    * of trail and error I have come up with a set filters to apply before
    * proceeding with the analysis.
    *
    * For a list of all of the things I have tried see:
    * ../scripts/00_Notebook_TOC.ipynb
    *
    * Here I implement the filtering strategy that I describe here:
    * ../scripts/ase_summary/ase_filters.ipynb
    *
    * Filters Include:
    *       Drop exonic regions that are always ambiguous in 100 genome simulation
    *       Drop exonic regions that have an APN < 25
    *       Drop exonic regions not present in at least 10% genotypes
    *       Drop exonic regions not in both environments
    *       Drop genotypes that show bias (median(q5_mean_theta) <= 0.4 or >= 0.6)
    *       Drop genotypes with too few exonic regions (<= 500)
    *
    * INPUT: CEGS.qsim_emp_theta_w_flag
    *
    * DATASET: CEGS.clean_ase_sbs     79,967 obs
    *          CEGS.clean_ase_stack   159,934 obs
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/ase_summary/exonic_regions_with_ai.csv
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_summarize_ase_filters.sas';

/* ASE Summary Stats */
    * A bunch of counts and summaries to get a better grasp on the ASE
    * results.
    *
    * 49 Genotypes remain after all filters
    * 1301 exonic regions have no AI in any genotype
    * 4090 exonic regions have AI in at least one genotype
    *
    * 14 exons have AI in more than 40 genotypes
    * 31 exons show AI in 100% of genotypes tested
    * 3004 exons show AI in 5 or fewer genotypes
    * 570 exons show AI in 10% or fewer genotypes tested
    *
    * After removing the 1303 exons that had no AI in any line and combining
    * mated and virgins..
    * w47 had the most exons with AI (n=1244, ~50% of exons);
    * r21 had the fewest exons with AI (n=143, ~24% of exons);
    * r857 had the highest percent of exons with AI (~56% exons tested);
    * r380 had the lowest percent of exons with AI (~13% exons tested);
    * Most genotypes appear to have ~25-30% of remaining exons with AI;
    *
    * There are 469 exonic regions where mated showed no AI and virgin had at
    * least 1 genotype with AI.
    *
    * There are 491 exonic regions where virgin showed no AI and mated had at
    * least 1 genotype with AI.
    *
    * There is only 1 exonic region where mated and virgin are more than 10
    * genotypes different. [F9836_SI, M=14, V=27]
    *
    * INPUT: CEGS.clean_ase_sbs
    *
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_summarize_ase_counts.sas';

/* ASE Cis-effects Enrichment Tests */
    * A set of enrichment tests to determine if there is an enrichment of
    * cis-effects.
    *
    * INPUT: CEGS.clean_ase_sbs
    *
    * DATASET: CEGS.cis_line_flags_for_go
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/ase_summary/cis_enrichment.rtf
    *
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_summarize_cis_enrichment.sas';
    
/*  Discordance between genotypes */
    * Want to see if there are groups of genotypes that are biased towards the
    * line and then others that are biased towards the tester.
    *
    * INPUT: CEGS.clean_ase_sbs
    *
    * DATASET: CEGS.discordant_genotypes
    *
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_summarize_ase_discordance_genotype.sas';

/* Discordance between mating status */
    * We are interested in exonic regions/genes that switch AI direction with
    * the environment. This could be something of biological interest. Pulls
    * out these exonic regions and does some basic counts.
    *
    * I don't think there is much here. Most fusions only show discordance
    * in 1 or 2 lines. Only [S43600_SI Pepck] and  [F22489_SI n-syb] had 6 and
    * 4 lines respectively. 5 of the 6 genotypes of [S43600_SI Pepck] were
    * mildly discordant (distance < 0.05). Only [F22489_SI n-syb] looked
    * potentially interesting, but with only 4 genotypes, I think it will be
    * hard to make any conclusions.
    *
    *       Are there any fusions with many genotypes showing discord?
    *       num_geno_with_discord    Frequency
    *       -----------------------------------
    *                      0              2508
    *                      1               170
    *                      2                31
    *                      3                 1
    *                      4                 1
    *                      6                 1
    *
    *
    *       Are there any genotypes with many fusions showing discord?
    *       num_fus_with_discord    Frequency
    *       ---------------------------------
    *                      0                8
    *                      1                8
    *                      2                8
    *                      3                4
    *                      4                1
    *                      5                3
    *                      6                6
    *                      7                1
    *                      8                1
    *                     10                1
    *                     12                3
    *                     14                2
    *                     21                2
    *                     23                1
    *
    * INPUT: CEGS.clean_ase_sbs
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/ase_summary/environment_discordance.csv
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_summarize_ase_discordance_mating_status.sas';

/* Sex Determination Summary */
    * Look only at the sex determination genes and determine what genes have
    * evidence for AI.
    *
    * The following sex det genes are not in the clean data:
    *       Rbp1 Rm62 Spf45 dsx her ix msl-2 tra tra2
    *
    * INPUT: CEGS.clean_ase_sbs
    *        GENELIST.sex_det
    *        DMEL551.si_fusion_2_gene_id
    *
    * DATASET: CEGS.sex_det_on_off
    *          CEGS.sex_det_geno_AI
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_summarize_ase_sex_det.sas';

/* Compare to Sergeys TF List */
    * Sergey has given me a list of genes that are affected by the eQTLs with
    * snps in TF binding sites. I want to merge to this list and see if any of
    * these have AI.
    *
    * INPUT:
    * !MCLAB/cegs_ase_paper/external_data/sergey_eqtl_transcription_factor_list.txt
    * DMEL551.fb551_si_fusion_2_gene_id 
    * 
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/ase_summarize_eqtl_tf.sas';
