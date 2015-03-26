/* 
 * REVISION:
 *      12/01/2011 Re-ran big merge because changed FDR flag for the contrast to .20
 *      12/27/2011 included all means.
 *      01/03/2012 changed the flybase 5.30 annotation a little, so need to adjust script for that.
 */

libname fru "!MCLAB/Fru_network/sasdata";
libname dmel530 "!MCLAB/useful_dmel_data/flybase530/sasdata";

/* prepare other datasets for merge */
proc sort data=dmel530.Fb530_si_fusions_unique_flagged;
    by fusion_id;
    run;

proc sort data=fru.model_fru_by_fusion;
    by fusion_id;
    run;

proc sort data=fru.all_meansv2;
     by fusion_id;
     run;

proc sort data=fru.flag_fdr_contrast_by_fusion;
    by fusion_id;
    run;

proc sort data=fru.flag_fail_normality_by_fusion;
    by fusion_id;
    run;

proc sort data=fru.flag_no_var;
    by fusion_id;
    run;

proc sort data=fru.flag_zero_mean;
    by fusion_id;
    run;

proc sort data=fru.flag_missing;
    by fusion_id;
    run;

data fru.results_by_fusion;
      merge dmel530.Fb530_si_fusions_unique_flagged (in=in1)
          fru.model_fru_by_fusion
          fru.all_meansv2
          fru.flag_fdr_contrast_by_fusion
          fru.flag_fail_normality_by_fusion
          fru.flag_no_var
          fru.flag_zero_mean
          fru.flag_missing;
    by fusion_id;
    if in1;
    run;

/* Make big results file with annotation */

proc sort data=fru.results_by_fusion;
    by fusion_id;
    run;

proc sort data=dmel530.fusions2go;
    by fusion_id;
    run;
    
proc contents data=fru.results_by_fusion varnum; run;

data fru.results_plus_gov2;
    retain Fusion_ID CHROM START END symbol_cat Exon_Gene_ID_cat Sequence_loc_cat exon_ID_cat Exon_Name_cat FBtrs_per_exon_cat FBgn_cat FBpp_cat FBtr_cat 
    max_fbtr_per_gene_symbol min_fbtr_per_gene_symbol_cat max_fbtr_per_gene_symbol_cat DF SS MS FValue ProbF AH_BerF_mean_logrpkm AH_BerF_sd_logrpkm 
    AH_BerM_mean_logrpkm AH_BerM_sd_logrpkm AH_CS_mean_logrpkm AH_CS_sd_logrpkm AH_CSFemale_mean_logrpkm AH_CSFemale_sd_logrpkm AH_Female_FruM_A__mean_logrpkm 
    AH_Female_FruM_A__sd_logrpkm AH_Female_FruM_B__mean_logrpkm AH_Female_FruM_B__sd_logrpkm AH_Female_FruM_C__mean_logrpkm AH_Female_FruM_C__sd_logrpkm 
    AH_FruP14_440_mean_logrpkm AH_FruP14_440_sd_logrpkm AH_FruW12_ChaM5_mean_logrpkm AH_FruW12_ChaM5_sd_logrpkm AH_Male_FruM_A__mean_logrpkm AH_Male_FruM_A__sd_logrpkm 
    AH_Male_FruM_B__mean_logrpkm AH_Male_FruM_B__sd_logrpkm AH_Male_FruM_C__mean_logrpkm AH_Male_FruM_C__sd_logrpkm AH_BerF_mean_rpkm AH_BerF_sd_rpkm AH_BerM_mean_rpkm 
    AH_BerM_sd_rpkm AH_CS_mean_rpkm AH_CS_sd_rpkm AH_CSFemale_mean_rpkm AH_CSFemale_sd_rpkm AH_Female_FruM_A__mean_rpkm AH_Female_FruM_A__sd_rpkm AH_Female_FruM_B__mean_rpkm 
    AH_Female_FruM_B__sd_rpkm AH_Female_FruM_C__mean_rpkm AH_Female_FruM_C__sd_rpkm AH_FruP14_440_mean_rpkm AH_FruP14_440_sd_rpkm AH_FruW12_ChaM5_mean_rpkm AH_FruW12_ChaM5_sd_rpkm 
    AH_Male_FruM_A__mean_rpkm AH_Male_FruM_A__sd_rpkm AH_Male_FruM_B__mean_rpkm AH_Male_FruM_B__sd_rpkm AH_Male_FruM_C__mean_rpkm AH_Male_FruM_C__sd_rpkm p_all_1 flag_p_all_1_05 
    p_all_2 flag_p_all_2_05 p_all_3 flag_p_all_3_05 p_all_4 flag_p_all_4_05 p_all_5 flag_p_all_5_05 p_all_6 flag_p_all_6_05 p_all_7 flag_p_all_7_05 p_all_8 flag_p_all_8_05 p_all_9 
    flag_p_all_9_05 p_all_10 flag_p_all_10_05 p_all_11 flag_p_all_11_05 p_all_12 flag_p_all_12_05 p_all_13 flag_p_all_13_05 p_all_14 flag_p_all_14_05 p_all_15 flag_p_all_15_05 
    p_all_16 flag_p_all_16_05 p_all_17 flag_p_all_17_05 p_all_18 flag_p_all_18_05 p_all_19 flag_p_all_19_05 p_all_20 flag_p_all_20_05 p_all_21 flag_p_all_21_05 p_all_22 
    flag_p_all_22_05 fdr_p_contrast_1 flag_fdr_p_contrast_1_20 fdr_p_contrast_2 flag_fdr_p_contrast_2_20 fdr_p_contrast_3 flag_fdr_p_contrast_3_20 fdr_p_contrast_4 
    flag_fdr_p_contrast_4_20 fdr_p_contrast_5 flag_fdr_p_contrast_5_20 fdr_p_contrast_6 flag_fdr_p_contrast_6_20 fdr_p_contrast_7 flag_fdr_p_contrast_7_20 
    fdr_p_contrast_8 flag_fdr_p_contrast_8_20 fdr_p_contrast_9 flag_fdr_p_contrast_9_20 fdr_p_contrast_10 flag_fdr_p_contrast_10_20 fdr_p_contrast_11 flag_fdr_p_contrast_11_20 
    fdr_p_contrast_12 flag_fdr_p_contrast_12_20 fdr_p_contrast_13 flag_fdr_p_contrast_13_20 fdr_p_contrast_14 flag_fdr_p_contrast_14_20 fdr_p_contrast_15 flag_fdr_p_contrast_15_20 
    fdr_p_contrast_16 flag_fdr_p_contrast_16_20 fdr_p_contrast_17 flag_fdr_p_contrast_17_20 fdr_p_contrast_18 flag_fdr_p_contrast_18_20 fdr_p_contrast_19 flag_fdr_p_contrast_19_20 
    fdr_p_contrast_20 flag_fdr_p_contrast_20_20 fdr_p_contrast_21 flag_fdr_p_contrast_21_20 fdr_p_contrast_22 flag_fdr_p_contrast_22_20 flag_fail_normality flag_no_var 
    flag_zero_mean flag_missing Genes_per_fusion Exons_per_fusion FBgns_per_fusion FBpps_per_fusion FBtrs_per_fusion Constitutive Common Alternative most_three_prime_exon 
    biopro_GO_ID_count_cat GO_number_biopro_cat_cat GO_biological_process_cat_cat molfunc_GO_ID_count_cat GO_number_molfunc_cat_cat GO_molecular_function_cat_cat 
    celcomp_GO_ID_count_cat GO_number_celcomp_cat_cat GO_cellular_component_cat_cat ;
    merge fru.results_by_fusion (in=in1) dmel530.fusions2go; 
    by fusion_id;
    if in1; 
    label df=; label ss=;
    label ms=; label fvalue=; 
    label probf=;
    run; 
