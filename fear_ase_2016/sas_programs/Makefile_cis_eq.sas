/*******************************************************************************
* Filename: Makefile_cis_eq.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Use the maren equations to come up with an estimate of
* cis-effects.
*******************************************************************************/

/* Libraries */
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);

/* Data Normalization */
    * For the maren equations, I am also going to drop exonic regions with less
    * than 10 genotypes. The maren equations make some assumptions about the
    * population level sums. Obvisouly the more genotypes that are present for
    * each fusions the better, but I am comfortable with as few as 10 genotypes.
    *
    * INPUT: CEGS.clean_ase_stack
    *
    * DATASET: CEGS.data_for_cis_eq
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/cis_eq_data_normalization.sas';

/* Maren Equations */
    * Sergey and Lauren developed a set of equtions found here:
    * 
    * Nuzhdin, S. V, Friesen, M. L., & McIntyre, L. M. (2012).
    * Genotype-phenotype mapping in a post-GWAS world. Trends in Genetics : TIG,
    * 28(9), 421–6. doi:10.1016/j.tig.2012.06.003
    * 
    * which potentially allow for the identificqtion of cis- and trans-effects.
    * Here I try using these qeustions and test if they give reasonable results.
    * 
    * Basics: For a given gene the expression level of Eii of allele i in F1
    * genotype i.
    * 
    * Eii=μ+Ci+(Ti+Tt)/2
    * 
    * Eti=μ+Ct+(Ti+Tt)/2
    * 
    * For each allele the cis- and trans-effects are deviations from the
    * population means, we expect that they will sum to zero:
    * 
    * ∑ni=1Ci=0
    * 
    * ∑ni=1Ti=0
    * 
    * Then the expected difference in expression between the Line and Tester
    * allele over the entire population is:
    * 
    * ∑ni=1Eti−Eiin
    * 
    * Which can be re-written as
    * 
    * ∑ni=1Ct−Cin=Ct
    * 
    * The cis-effect of allele i can be estimated by:
    * 
    * ˆCi=ˆEii−ˆEti+ˆCt
    *
    * INPUT: CEGS.data_for_cis_eq
    *
    * DATASET: CEGS.cis_eq
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/cis_effects/cis_eq_full_output_w_ai_calls.csv
    *
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/cis_eq.sas';

/* Export data matrices for clustering */
    *
    * INPUT: CEGS.cis_eq
    *
    * DATASET: CEGS.cis_line_effect_mated
    *          CEGS.cis_line_effect_virgin
    *          CEGS.trans_line_effect_mated
    *          CEGS.trans_line_effect_virgin
    *
    * FILE: !MCLAB/cegs_ase_paper/pipeline_output/cis_effects/cis_line_effect_mated.csv
    *       !MCLAB/cegs_ase_paper/pipeline_output/cis_effects/cis_line_effect_virgin.csv
    *       !MCLAB/cegs_ase_paper/pipeline_output/cis_effects/trans_line_effect_mated.csv
    *       !MCLAB/cegs_ase_paper/pipeline_output/cis_effects/trans_line_effect_virgin.csv
    *       !MCLAB/cegs_ase_paper/pipeline_output/cis_effects/mmc/cis_line_effect_for_mmc_mated.csv
    *       !MCLAB/cegs_ase_paper/pipeline_output/cis_effects/mmc/cis_line_effect_for_mmc_virgin.csv
    ;
    %include '!MCLAB/cegs_ase_paper/sas_programs/cis_eq_for_clustering.sas';
