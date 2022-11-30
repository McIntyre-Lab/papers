#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=(
        "Prepare datasets for scatter plots of female/male ratios from exon fragments significant "
        "for DE in D. melanogaster vs. D. simulans."
        )
    )

    # Input data
    parser.add_argument(
        "-mf",
        "--mel-feature",
        dest="melFeature",
        required=True,
        help="D. melanogaster CSV with exon fragment f/m ratios and DE significance flags."
    )
    parser.add_argument(
        "-ma",
        "--mel-annotation",
        dest="melAnnot",
        required=True,
        help="D. melanogaster annotation of exon fragment to gene_id (FBgn)."
    )
    parser.add_argument(
        "-sf",
        "--sim-feature",
        dest="simFeature",
        required=True,
        help="D. simulans CSV with exon fragment f/m ratios and DE significance flags."
    )
    parser.add_argument(
        "-sa",
        "--sim-annotation",
        dest="simAnnot",
        required=True,
        help="D. simulans annotation of exon fragment to gene_id (FBgn)."
    )
    parser.add_argument(
        "-o",
        "--ortho",
        dest="ortho",
        required=True,
        help="CSV of D. melanogaster to D. simulans gene-level ortholog results."
    )

    # Output data
    parser.add_argument(
        "-d",
        "--out-directory",
        dest="outDir",
        required=True,
        help="Output directory for CSV files."
    )

    args = parser.parse_args()
    return args

def main():
    # Get input files
    melFeatureDf = pd.read_csv(args.melFeature, low_memory=False)
    melAnnotDf = pd.read_csv(args.melAnnot, low_memory=False)
    simFeatureDf = pd.read_csv(args.simFeature, low_memory=False)
    simAnnotDf = pd.read_csv(args.simAnnot, low_memory=False)
    orthoDf = pd.read_csv(args.ortho, low_memory=False)
    
    # melFeatureDf = pd.read_csv("mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq/model_output/mel_frag_flags_kitchen_sink.csv", low_memory=False)
    # melAnnotDf = pd.read_csv("mclab/SHARE/McIntyre_Lab/useful_dmel_data/flybase617/event_analysis_annotations/150bp_annotations/dmel617_exon_fragment_annotations.csv", low_memory=False)
    # simFeatureDf = pd.read_csv("mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq/model_output/sim_frag_flags_kitchen_sink.csv", low_memory=False)
    # simAnnotDf = simAnnotDf = pd.read_csv("mclab/SHARE/McIntyre_Lab/useful_dsim_data/flybase202/event_analysis_annotations/dsim202_annotations_150bp_reads/dsim202_150bp_exon_fragment_annotations.csv", low_memory=False)
    # orthoDf = pd.read_csv("/Users/adalena/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/supp_files/dmel_dsim_ortholog_chip_rna_flags.csv", low_memory=False)

    # Select only significant features for each species
    melSigFeatureDf = melFeatureDf[
        melFeatureDf["flag_ttest_1_trend_F"] + melFeatureDf["flag_ttest_1_trend_M"] > 0
    ].copy()
    simSigFeatureDf = simFeatureDf[
        simFeatureDf["flag_ttest_1_trend_F"] + simFeatureDf["flag_ttest_1_trend_M"] > 0
    ].copy()

    # Make M/F ratio (F/M ratio already present)
    melSigFeatureDf["ratio_avg_m_f_apn_uq_ff"] = (
        melSigFeatureDf["avg_m_apn_uq_ff"] / melSigFeatureDf["avg_f_apn_uq_ff"]
    )
    simSigFeatureDf["ratio_avg_m_f_apn_uq_ff"] = (
        simSigFeatureDf["avg_m_apn_uq_ff"] / simSigFeatureDf["avg_f_apn_uq_ff"]
    )

    # Select only one-to-one orthologs
    one2oneDf = orthoDf[orthoDf["flag_one2one_ortholog"]==1].copy()
    del(orthoDf)

    # Merge fragment annotation into feature flags and ratios to get only fragments
    melSigFeatureGeneDf = pd.merge(
        melSigFeatureDf,
        melAnnotDf[["fragment_id", "gene_id"]],
        how = "outer",
        left_on = "featureID",
        right_on = "fragment_id",
        indicator = "merge_check",
        validate = "1:1"
    )

    simSigFeatureGeneDf = pd.merge(
        simSigFeatureDf,
        simAnnotDf[["fragment_id", "gene_id"]],
        how = "outer",
        left_on = "featureID",
        right_on = "fragment_id",
        indicator = "merge_check",
        validate = "1:1"
    )

    # Select only fragments and drop unneccessary columns
    melSigFragGeneDf = melSigFeatureGeneDf[
            melSigFeatureGeneDf["merge_check"]=="both"
        ][[
            "gene_id",
            "featureID",
            "avg_f_apn_uq_ff",
            "avg_m_apn_uq_ff",
            "ratio_avg_f_m_apn_uq_ff",
            "ratio_avg_m_f_apn_uq_ff",
            "flag_ttest_1_trend_F",
            "flag_ttest_1_trend_M"]].copy()
    
    simSigFragGeneDf = simSigFeatureGeneDf[
            simSigFeatureGeneDf["merge_check"]=="both"
        ][[
            "gene_id",
            "featureID",
            "avg_f_apn_uq_ff",
            "avg_m_apn_uq_ff",
            "ratio_avg_f_m_apn_uq_ff",
            "ratio_avg_m_f_apn_uq_ff",
            "flag_ttest_1_trend_F",
            "flag_ttest_1_trend_M"]].copy()

    # Check that all other features were 3UTR, 5UTR, intron, or TSS
    checkMel = len(
        melSigFeatureGeneDf[
                (melSigFeatureGeneDf["merge_check"]=="left_only")
                & (~melSigFeatureGeneDf["featureID"].str.contains("3UTR", na=False))
                & (~melSigFeatureGeneDf["featureID"].str.contains("5UTR", na=False))
                & (~melSigFeatureGeneDf["featureID"].str.contains("intron", na=False))
                & (~melSigFeatureGeneDf["featureID"].str.contains("TSS", na=False))
            ].copy().drop(columns=["merge_check"])
    )
    if checkMel != 0:
        print("WARNING: Unexpected features detected in D. melanogaster featureID column.")
    del(melFeatureDf, melSigFeatureDf, melSigFeatureGeneDf)

    checkSim = len(
        simSigFeatureGeneDf[
                (simSigFeatureGeneDf["merge_check"]=="left_only")
                & (~simSigFeatureGeneDf["featureID"].str.contains("3UTR", na=False))
                & (~simSigFeatureGeneDf["featureID"].str.contains("5UTR", na=False))
                & (~simSigFeatureGeneDf["featureID"].str.contains("intron", na=False))
                & (~simSigFeatureGeneDf["featureID"].str.contains("TSS", na=False))
            ].copy().drop(columns=["merge_check"])
    )
    if checkSim != 0:
        print("WARNING: Unexpected features detected in D. simulans featureID column.")
    del(simFeatureDf, simSigFeatureDf, simSigFeatureGeneDf)

    # Summarize gene using the min and max ratio of the fragments within the gene
    melSigGeneDf = melSigFragGeneDf.groupby("gene_id").agg({
        "ratio_avg_f_m_apn_uq_ff": ["min", "max"],
        "ratio_avg_m_f_apn_uq_ff": ["min", "max"]}
    )
    melSigGeneCols = [c[1]+"_"+c[0] for c in melSigGeneDf.columns]
    melSigGeneDf.columns = melSigGeneCols
    melSigGeneDf = melSigGeneDf.reset_index()

    simSigGeneDf = simSigFragGeneDf.groupby("gene_id").agg({
        "ratio_avg_f_m_apn_uq_ff": ["min", "max"],
        "ratio_avg_m_f_apn_uq_ff": ["min", "max"]}
    )
    simSigGeneCols = [c[1]+"_"+c[0] for c in simSigGeneDf.columns]
    simSigGeneDf.columns = simSigGeneCols
    simSigGeneDf = simSigGeneDf.reset_index()

    # Get expression flags from ortholog gene-level results file
    expFlagLst = ['flag_U_mel', 'flag_U_sim', 'flag_M_mel', 'flag_M_sim', 'flag_M_mel_only',
                  'flag_M_sim_only', 'flag_M_mel_U_sim', 'flag_M_sim_U_mel', 'flag_MP_mel',
                  'flag_MP_sim', 'flag_F_mel', 'flag_F_sim', 'flag_F_mel_only', 'flag_F_sim_only',
                  'flag_F_mel_U_sim', 'flag_F_sim_U_mel', 'flag_M_mel_F_sim', 'flag_M_sim_F_mel',
                  'flag_FP_mel', 'flag_FP_sim', 'flag_MF_mel', 'flag_MF_sim', 'flag_MF_mel_only',
                  'flag_MF_sim_only', 'flag_MF_mel_U_sim', 'flag_MF_sim_U_mel', 'flag_MFP_mel',
                  'flag_MFP_sim', 'flag_sex_bias_both_species', 'flag_conserved_sex_bias',
                  'flag_conserved_male_bias', 'flag_conserved_female_bias',
                  'flag_conserved_male_and_female', 'flag_conserved_unbiased', 
                  'flag_sex_bias_mel_only', 'flag_sex_bias_sim_only', 'flag_switch_M_U',
                  'flag_switch_F_U', 'flag_switch_M_F', 'flag_diverged_sex_bias',
                  'flag_has_male_k4_mel', 'flag_has_male_k4_sim', 'flag_has_male_k4_both_species',
                  'flag_has_female_k4_mel', 'flag_has_female_k4_sim', 'flag_has_female_k4_both_species',
                  'flag_has_male_k27_mel', 'flag_has_male_k27_sim', 'flag_has_male_k27_both_species',
                  'flag_has_female_k27_mel', 'flag_has_female_k27_sim', 'flag_has_female_k27_both_species',]
    GOflagLst = ['biological_process', 'goterm_biol_process', 'molecular_function',
                 'goterm_mol_function', 'cellular_component', 'goterm_cell_component',]
    mappingFlagLst = ['flag_map_better_2_mel', 'flag_map_better_2_sim']


    # Merge mel with orthologous genes
    melOrthoDf = pd.merge(
        melSigGeneDf,
        one2oneDf[["mel_geneID", "mel_geneSymbol", "sim_geneID", "sim_geneSymbol", "xsome"] + expFlagLst + GOflagLst + mappingFlagLst],
        how = "outer",
        left_on = "gene_id",
        right_on = "mel_geneID",
        indicator = "merge_check",
        validate = "1:1",
    )
    melOrthoDf = melOrthoDf.rename(columns={
        "gene_id": "mel_gene_id",
        "featureID": "mel_featureID",
        "avg_f_apn_uq_ff": "mel_avg_f_apn_uq_ff",
        "avg_m_apn_uq_ff": "mel_avg_m_apn_uq_ff",
        "min_ratio_avg_f_m_apn_uq_ff": "mel_min_ratio_avg_f_m_apn_uq_ff",
        "max_ratio_avg_f_m_apn_uq_ff": "mel_max_ratio_avg_f_m_apn_uq_ff",
        "min_ratio_avg_m_f_apn_uq_ff": "mel_min_ratio_avg_m_f_apn_uq_ff",
        "max_ratio_avg_m_f_apn_uq_ff": "mel_max_ratio_avg_m_f_apn_uq_ff",
        "flag_ttest_1_trend_F": "mel_flag_ttest_1_trend_F",
        "flag_ttest_1_trend_M": "mel_flag_ttest_1_trend_M",
    })
    # Drop genes without an ortholog or significant fragment (keep "both")
    melOrthoDf = melOrthoDf[
            melOrthoDf["merge_check"]=="both"
        ].drop(columns=["merge_check", "mel_geneID"]).copy()

    # Merge sim with mel and orthologs
    melSimOrthoDf = pd.merge(
        simSigGeneDf,
        melOrthoDf,
        how = "outer",
        left_on = "gene_id",
        right_on = "sim_geneID",
        indicator = "merge_check",
        validate = "1:1"
    )
    melSimOrthoDf = melSimOrthoDf.rename(columns={
        "gene_id": "sim_gene_id",
        "featureID": "sim_featureID",
        "avg_f_apn_uq_ff": "sim_avg_f_apn_uq_ff",
        "avg_m_apn_uq_ff": "sim_avg_m_apn_uq_ff",
        "min_ratio_avg_f_m_apn_uq_ff": "sim_min_ratio_avg_f_m_apn_uq_ff",
        "max_ratio_avg_f_m_apn_uq_ff": "sim_max_ratio_avg_f_m_apn_uq_ff",
        "min_ratio_avg_m_f_apn_uq_ff": "sim_min_ratio_avg_m_f_apn_uq_ff",
        "max_ratio_avg_m_f_apn_uq_ff": "sim_max_ratio_avg_m_f_apn_uq_ff",
        "flag_ttest_1_trend_F": "sim_flag_ttest_1_trend_F",
        "flag_ttest_1_trend_M": "sim_flag_ttest_1_trend_M",
    })
    # Drop genes without an ortholog or significant fragment (keep "both")
    melSimOrthoDf = melSimOrthoDf[
            melSimOrthoDf["merge_check"]=="both"
        ].drop(columns=["merge_check", "sim_geneID"]).copy()

    # For genes that are male and female-biased expression in both species
    #   remove comparison of fragments biased towards the other sex
    # melSimOrthoDf = melSimOrthoDf[
    #     (melSimOrthoDf["flag_MF_mel"] + melSimOrthoDf["flag_MF_sim"] != 2)
    #     | (
    #         (melSimOrthoDf["flag_MF_mel"] + melSimOrthoDf["flag_MF_sim"] == 2)
    #         & (melSimOrthoDf["mel_flag_ttest_1_trend_F"] == melSimOrthoDf["sim_flag_ttest_1_trend_F"])
    #         & (melSimOrthoDf["mel_flag_ttest_1_trend_M"] == melSimOrthoDf["sim_flag_ttest_1_trend_M"])
    #     )
    # ].copy()

    # Drop fragments where the male or female average was 0 in either species (set to 10^-6)
    # melSimOrthoDf = melSimOrthoDf[
    #     (melSimOrthoDf["mel_avg_f_apn_uq_ff"] != 0.000001)
    #     & (melSimOrthoDf["mel_avg_m_apn_uq_ff"] != 0.000001)
    #     & (melSimOrthoDf["sim_avg_m_apn_uq_ff"] != 0.000001)
    #     & (melSimOrthoDf["sim_avg_f_apn_uq_ff"] != 0.000001)
    # ].copy()

    # Get conserved male-biased and female-biased genes
    MconservDf = melSimOrthoDf[
        (melSimOrthoDf["flag_conserved_male_bias"]==1)
        & (melSimOrthoDf["xsome"].isin(["X", "A"]))
    ].copy()
    FconservDf = melSimOrthoDf[
        (melSimOrthoDf["flag_conserved_female_bias"]==1)
        & (melSimOrthoDf["xsome"].isin(["X", "A"]))
    ].copy()

    # Make ratios:
    # both constrained between 0 and 1 w/ lowest bias at 0 (1 - min M/F for female-biased, 1 - min F/M for male-biased)
    MconservDf["one_minus_F_M_ratio_mel"] = 1 - MconservDf["mel_min_ratio_avg_f_m_apn_uq_ff"]
    MconservDf["one_minus_F_M_ratio_sim"] = 1 - MconservDf["sim_min_ratio_avg_f_m_apn_uq_ff"]
    FconservDf["one_minus_M_F_ratio_mel"] = 1 - FconservDf["mel_min_ratio_avg_m_f_apn_uq_ff"]
    FconservDf["one_minus_M_F_ratio_sim"] = 1 - FconservDf["sim_min_ratio_avg_m_f_apn_uq_ff"]

    # Output convserved ratio CSV
    chromFlagLst = ['flag_has_male_k4_mel', 'flag_has_male_k4_sim', 'flag_has_female_k4_mel',
                    'flag_has_female_k4_sim', 'flag_has_male_k27_mel', 'flag_has_male_k27_sim',
                    'flag_has_female_k27_mel', 'flag_has_female_k27_sim']
    MconservDf[[
        "mel_gene_id",
        "mel_geneSymbol",
        "sim_gene_id",
        "sim_geneSymbol",
        "one_minus_F_M_ratio_mel",
        "one_minus_F_M_ratio_sim"] + GOflagLst + chromFlagLst + mappingFlagLst].to_csv(args.outDir + "/ortho_MF_ratio_conserved_M.csv", index=False)
    FconservDf[[
        "mel_gene_id",
        "mel_geneSymbol",
        "sim_gene_id",
        "sim_geneSymbol",
        "one_minus_M_F_ratio_mel",
        "one_minus_M_F_ratio_sim"] + GOflagLst + chromFlagLst + mappingFlagLst].to_csv(args.outDir + "/ortho_MF_ratio_conserved_F.csv", index=False)


    # Get divergent sex biased genes
    MmelFsimDf = melSimOrthoDf[
        (melSimOrthoDf["flag_M_mel_F_sim"]==1)
        & (melSimOrthoDf["xsome"].isin(["X", "A"]))
    ].copy()
    MmelFsimDf["one_minus_F_M_ratio_mel"] = 1 - MmelFsimDf["mel_min_ratio_avg_f_m_apn_uq_ff"]
    MmelFsimDf["one_minus_M_F_ratio_sim"] = 1 - MmelFsimDf["sim_min_ratio_avg_m_f_apn_uq_ff"]

    FmelMsimDf = melSimOrthoDf[
        (melSimOrthoDf["flag_M_sim_F_mel"]==1)
        & (melSimOrthoDf["xsome"].isin(["X", "A"]))
    ].copy()
    FmelMsimDf["one_minus_M_F_ratio_mel"] = 1 - FmelMsimDf["mel_min_ratio_avg_m_f_apn_uq_ff"]
    FmelMsimDf["one_minus_F_M_ratio_sim"] = 1 - FmelMsimDf["sim_min_ratio_avg_f_m_apn_uq_ff"]

    # Output divergent sex biased genes
    MmelFsimDf[[
        "mel_gene_id",
        "mel_geneSymbol",
        "sim_gene_id",
        "sim_geneSymbol",
        "one_minus_F_M_ratio_mel",
        "one_minus_M_F_ratio_sim"] + GOflagLst + chromFlagLst + mappingFlagLst].to_csv(args.outDir + "/ortho_MF_ratio_M_mel_F_sim.csv", index=False)
    FmelMsimDf[[
        "mel_gene_id",
        "mel_geneSymbol",
        "sim_gene_id",
        "sim_geneSymbol",
        "one_minus_M_F_ratio_mel",
        "one_minus_F_M_ratio_sim"] + GOflagLst + chromFlagLst + mappingFlagLst].to_csv(args.outDir + "/ortho_MF_ratio_F_mel_M_sim.csv", index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
