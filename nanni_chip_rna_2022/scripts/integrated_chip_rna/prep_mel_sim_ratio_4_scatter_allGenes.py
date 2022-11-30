#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=(
        "Prepare datasets for scatter plots of female/male ratios from exon fragments "
        "either unbiased or sex-biased in D. melanogaster vs. D. simulans."
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


    # Select only one-to-one orthologs
    one2oneDf = orthoDf[orthoDf["flag_one2one_ortholog"]==1].copy()
    del(orthoDf)

    # Make M/F ratio (F/M ratio already present)
    melFeatureDf["ratio_avg_m_f_apn_uq_ff"] = (
        melFeatureDf["avg_m_apn_uq_ff"] / melFeatureDf["avg_f_apn_uq_ff"]
    )
    simFeatureDf["ratio_avg_m_f_apn_uq_ff"] = (
        simFeatureDf["avg_m_apn_uq_ff"] / simFeatureDf["avg_f_apn_uq_ff"]
    )

    # Get ratios for signifiant features only
    melFeatureDf["sig_ratio_avg_f_m_apn_uq_ff"] = np.where(
        melFeatureDf["flag_ttest_1_trend_F"] + melFeatureDf["flag_ttest_1_trend_M"] > 0,
        melFeatureDf["ratio_avg_f_m_apn_uq_ff"],
        np.nan
    )
    simFeatureDf["sig_ratio_avg_f_m_apn_uq_ff"] = np.where(
        simFeatureDf["flag_ttest_1_trend_F"] + simFeatureDf["flag_ttest_1_trend_M"] > 0,
        simFeatureDf["ratio_avg_f_m_apn_uq_ff"],
        np.nan
    )
    melFeatureDf["sig_ratio_avg_m_f_apn_uq_ff"] = np.where(
        melFeatureDf["flag_ttest_1_trend_F"] + melFeatureDf["flag_ttest_1_trend_M"] > 0,
        melFeatureDf["ratio_avg_m_f_apn_uq_ff"],
        np.nan
    )
    simFeatureDf["sig_ratio_avg_m_f_apn_uq_ff"] = np.where(
        simFeatureDf["flag_ttest_1_trend_F"] + simFeatureDf["flag_ttest_1_trend_M"] > 0,
        simFeatureDf["ratio_avg_m_f_apn_uq_ff"],
        np.nan
    )

    # Merge fragment annotation into feature flags and ratios to get only fragments
    melFeatureGeneDf = pd.merge(
        melFeatureDf,
        melAnnotDf[["fragment_id", "gene_id"]],
        how = "outer",
        left_on = "featureID",
        right_on = "fragment_id",
        indicator = "merge_check",
        validate = "1:1"
    )

    simFeatureGeneDf = pd.merge(
        simFeatureDf,
        simAnnotDf[["fragment_id", "gene_id"]],
        how = "outer",
        left_on = "featureID",
        right_on = "fragment_id",
        indicator = "merge_check",
        validate = "1:1"
    )

    # Select only fragments and drop unneccessary columns
    melFragGeneDf = melFeatureGeneDf[
            melFeatureGeneDf["merge_check"]=="both"
        ][[
            "gene_id",
            "featureID",
            "avg_f_apn_uq_ff",
            "avg_m_apn_uq_ff",
            "ratio_avg_f_m_apn_uq_ff",
            "ratio_avg_m_f_apn_uq_ff",
            "sig_ratio_avg_f_m_apn_uq_ff",
            "sig_ratio_avg_m_f_apn_uq_ff",
            "flag_ttest_1_trend_F",
            "flag_ttest_1_trend_M"]].copy()
    
    simFragGeneDf = simFeatureGeneDf[
            simFeatureGeneDf["merge_check"]=="both"
        ][[
            "gene_id",
            "featureID",
            "avg_f_apn_uq_ff",
            "avg_m_apn_uq_ff",
            "ratio_avg_f_m_apn_uq_ff",
            "ratio_avg_m_f_apn_uq_ff",
            "sig_ratio_avg_f_m_apn_uq_ff",
            "sig_ratio_avg_m_f_apn_uq_ff",
            "flag_ttest_1_trend_F",
            "flag_ttest_1_trend_M"]].copy()

    # Check that all other features were 3UTR, 5UTR, intron, or TSS
    checkMel = len(
        melFeatureGeneDf[
                (melFeatureGeneDf["merge_check"]=="left_only")
                & (~melFeatureGeneDf["featureID"].str.contains("3UTR", na=False))
                & (~melFeatureGeneDf["featureID"].str.contains("5UTR", na=False))
                & (~melFeatureGeneDf["featureID"].str.contains("intron", na=False))
                & (~melFeatureGeneDf["featureID"].str.contains("TSS", na=False))
            ].copy().drop(columns=["merge_check"])
    )
    if checkMel != 0:
        print("WARNING: Unexpected features detected in D. melanogaster featureID column.")
    del(melFeatureDf, melFeatureGeneDf)

    checkSim = len(
        simFeatureGeneDf[
                (simFeatureGeneDf["merge_check"]=="left_only")
                & (~simFeatureGeneDf["featureID"].str.contains("3UTR", na=False))
                & (~simFeatureGeneDf["featureID"].str.contains("5UTR", na=False))
                & (~simFeatureGeneDf["featureID"].str.contains("intron", na=False))
                & (~simFeatureGeneDf["featureID"].str.contains("TSS", na=False))
            ].copy().drop(columns=["merge_check"])
    )
    if checkSim != 0:
        print("WARNING: Unexpected features detected in D. simulans featureID column.")
    del(simFeatureDf, simFeatureGeneDf)


    # Summarize gene using the min, max, and mean ratio of the fragments within the gene
    melSigGeneDf = melFragGeneDf.groupby("gene_id").agg({
        "ratio_avg_f_m_apn_uq_ff": ["min", "max", "mean"],
        "ratio_avg_m_f_apn_uq_ff": ["min", "max", "mean"],
        "sig_ratio_avg_f_m_apn_uq_ff": ["min", "max", "mean"],
        "sig_ratio_avg_m_f_apn_uq_ff": ["min", "max", "mean"]}
    )
    melSigGeneCols = [c[1]+"_"+c[0] for c in melSigGeneDf.columns]
    melSigGeneDf.columns = melSigGeneCols
    melSigGeneDf = melSigGeneDf.reset_index()

    simSigGeneDf = simFragGeneDf.groupby("gene_id").agg({
        "ratio_avg_f_m_apn_uq_ff": ["min", "max", "mean"],
        "ratio_avg_m_f_apn_uq_ff": ["min", "max", "mean"],
        "sig_ratio_avg_f_m_apn_uq_ff": ["min", "max", "mean"],
        "sig_ratio_avg_m_f_apn_uq_ff": ["min", "max", "mean"]}
    )
    simSigGeneCols = [c[1]+"_"+c[0] for c in simSigGeneDf.columns]
    simSigGeneDf.columns = simSigGeneCols
    simSigGeneDf = simSigGeneDf.reset_index()


    # Get expression flags from ortholog gene-level results file
    expFlagLst = ['flag_expressed_mel','flag_expressed_sim', 'flag_sex_limited_mel',
                  'flag_sex_limited_sim', 'flag_U_mel', 'flag_U_sim', 'flag_M_mel',
                  'flag_M_sim', 'flag_M_mel_only', 'flag_M_sim_only', 'flag_M_mel_U_sim',
                  'flag_M_sim_U_mel', 'flag_MP_mel', 'flag_MP_sim', 'flag_F_mel',
                  'flag_F_sim', 'flag_F_mel_only', 'flag_F_sim_only', 'flag_F_mel_U_sim',
                  'flag_F_sim_U_mel', 'flag_M_mel_F_sim', 'flag_M_sim_F_mel',
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
        one2oneDf[["mel_geneID", "mel_geneSymbol", "sim_geneID", "sim_geneSymbol", "xsome", "sim_xsome"] + expFlagLst + GOflagLst + mappingFlagLst],
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
        "mean_ratio_avg_f_m_apn_uq_ff": "mel_mean_ratio_avg_f_m_apn_uq_ff",
        "min_ratio_avg_m_f_apn_uq_ff": "mel_min_ratio_avg_m_f_apn_uq_ff",
        "max_ratio_avg_m_f_apn_uq_ff": "mel_max_ratio_avg_m_f_apn_uq_ff",
        "mean_ratio_avg_m_f_apn_uq_ff": "mel_mean_ratio_avg_m_f_apn_uq_ff",
        "min_sig_ratio_avg_f_m_apn_uq_ff": "mel_min_sig_ratio_avg_f_m_apn_uq_ff",
        "max_sig_ratio_avg_f_m_apn_uq_ff": "mel_max_sig_ratio_avg_f_m_apn_uq_ff",
        "mean_sig_ratio_avg_f_m_apn_uq_ff": "mel_mean_sig_ratio_avg_f_m_apn_uq_ff",
        "min_sig_ratio_avg_m_f_apn_uq_ff": "mel_min_sig_ratio_avg_m_f_apn_uq_ff",
        "max_sig_ratio_avg_m_f_apn_uq_ff": "mel_max_sig_ratio_avg_m_f_apn_uq_ff",
        "mean_sig_ratio_avg_m_f_apn_uq_ff": "mel_mean_sig_ratio_avg_m_f_apn_uq_ff",
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
        "mean_ratio_avg_f_m_apn_uq_ff": "sim_mean_ratio_avg_f_m_apn_uq_ff",
        "min_ratio_avg_m_f_apn_uq_ff": "sim_min_ratio_avg_m_f_apn_uq_ff",
        "max_ratio_avg_m_f_apn_uq_ff": "sim_max_ratio_avg_m_f_apn_uq_ff",
        "mean_ratio_avg_m_f_apn_uq_ff": "sim_mean_ratio_avg_m_f_apn_uq_ff",
        "min_sig_ratio_avg_f_m_apn_uq_ff": "sim_min_sig_ratio_avg_f_m_apn_uq_ff",
        "max_sig_ratio_avg_f_m_apn_uq_ff": "sim_max_sig_ratio_avg_f_m_apn_uq_ff",
        "mean_sig_ratio_avg_f_m_apn_uq_ff": "sim_mean_sig_ratio_avg_f_m_apn_uq_ff",
        "min_sig_ratio_avg_m_f_apn_uq_ff": "sim_min_sig_ratio_avg_m_f_apn_uq_ff",
        "max_sig_ratio_avg_m_f_apn_uq_ff": "sim_max_sig_ratio_avg_m_f_apn_uq_ff",
        "mean_sig_ratio_avg_m_f_apn_uq_ff": "sim_mean_sig_ratio_avg_m_f_apn_uq_ff",
        "flag_ttest_1_trend_F": "sim_flag_ttest_1_trend_F",
        "flag_ttest_1_trend_M": "sim_flag_ttest_1_trend_M",
    })
    # Drop genes without an ortholog
    melSimOrthoDf = melSimOrthoDf[
            melSimOrthoDf["merge_check"]=="both"
        ].drop(columns=["merge_check", "sim_geneID"]).copy()

    # Drop genes not on X or autosomes (flip X/A across the species) or not expression in both species
    melSimOrthoDf = melSimOrthoDf[
        (melSimOrthoDf["xsome"].isin(["X", "A"]))
        & (melSimOrthoDf["xsome"] == melSimOrthoDf["sim_xsome"])
        & (melSimOrthoDf["flag_expressed_mel"]+melSimOrthoDf["flag_expressed_sim"]==2)
    ]

    # Drop genes that are sex-limited in either species
    print(
        "{} genes are sex-limited in a least one species ({} in mel, {} in sim) "
        "- removed from plotting".format(
            len(melSimOrthoDf[melSimOrthoDf["flag_sex_limited_mel"]+melSimOrthoDf["flag_sex_limited_sim"]>0]),
            len(melSimOrthoDf[melSimOrthoDf["flag_sex_limited_mel"]>0]),
            len(melSimOrthoDf[melSimOrthoDf["flag_sex_limited_sim"]>0]),            
        )
    )
    melSimOrthoDf = melSimOrthoDf[
        melSimOrthoDf["flag_sex_limited_mel"]+melSimOrthoDf["flag_sex_limited_sim"]==0
    ]

    # Drop genes that are male and female expression in at least one of the species
    print(
        "{} genes classified as male- and female-biased ({} in mel, {} in sim) "
        "- removed from plotting".format(
            len(melSimOrthoDf[melSimOrthoDf["flag_MF_mel"]+melSimOrthoDf["flag_MF_sim"]>0]),
            len(melSimOrthoDf[melSimOrthoDf["flag_MF_mel"]>0]),
            len(melSimOrthoDf[melSimOrthoDf["flag_MF_sim"]>0])
        )
    )
    melSimOrthoDf = melSimOrthoDf[
        melSimOrthoDf["flag_MF_mel"] + melSimOrthoDf["flag_MF_sim"] == 0
    ].copy()

    # Count gene categories that will be plotted
    print("{} total genes:\n".format(len(melSimOrthoDf))
        + "{} conserved male-biased genes".format(int(melSimOrthoDf["flag_conserved_male_bias"].sum()))
        + "\n{} conserved female-biased genes".format(int(melSimOrthoDf["flag_conserved_female_bias"].sum()))
        + "\n{} conserved unbiased genes".format(int(melSimOrthoDf["flag_conserved_unbiased"].sum()))
        + "\n{} male-biased in mel and unbiased in sim".format(int(melSimOrthoDf["flag_M_mel_U_sim"].sum()))
        + "\n{} male-biased in sim and unbiased in mel".format(int(melSimOrthoDf["flag_M_sim_U_mel"].sum()))
        + "\n{} female-biased in mel and unbiased in sim".format(int(melSimOrthoDf["flag_F_mel_U_sim"].sum()))
        + "\n{} female-biased in sim and unbiased in mel".format(int(melSimOrthoDf["flag_F_sim_U_mel"].sum()))
        + "\n{} male-biased in mel and female-biased in sim".format(int(melSimOrthoDf["flag_M_mel_F_sim"].sum()))
        + "\n{} male-biased in sim and female-biased in mel".format(int(melSimOrthoDf["flag_M_sim_F_mel"].sum()))
    )
    # Total 9294

    # Set up ratios for each group of genes
    # Top right quadrant will be conserved male (1 - F/M, both species)
    # Bottom left quadrant will be conserved female (M/F - 1, both species)
    # Bottom right quadrant will be mel male and simulans female
    # Bottom left quadrant will be sim male and mel female

    # For unbiased gense where the mean ratios of M/F and F/M are both > 1:
    #   plot on the "high end" of the sex with more, e.g.,
    #   if M/F is higher than F/M set ratio to be 1-0 = 1
    #   if F/M is higher than M/F set ratio to be 0-1 = -1

    # Mel ratios
    melRatioConditions = [
        melSimOrthoDf["flag_M_mel"] == 1,
        melSimOrthoDf["flag_F_mel"] == 1,
        (melSimOrthoDf["flag_U_mel"] == 1)
            & (melSimOrthoDf["mel_mean_ratio_avg_f_m_apn_uq_ff"] < 1)
            & (melSimOrthoDf["mel_mean_ratio_avg_m_f_apn_uq_ff"] > 1),
        (melSimOrthoDf["flag_U_mel"] == 1)
            & (melSimOrthoDf["mel_mean_ratio_avg_m_f_apn_uq_ff"] < 1)
            & (melSimOrthoDf["mel_mean_ratio_avg_f_m_apn_uq_ff"] > 1),
        (melSimOrthoDf["flag_U_mel"] == 1)
            & (melSimOrthoDf["mel_mean_ratio_avg_m_f_apn_uq_ff"] > melSimOrthoDf["mel_mean_ratio_avg_f_m_apn_uq_ff"]),
        (melSimOrthoDf["flag_U_mel"] == 1)
            & (melSimOrthoDf["mel_mean_ratio_avg_m_f_apn_uq_ff"] < melSimOrthoDf["mel_mean_ratio_avg_f_m_apn_uq_ff"])
    ]
    melRatioChoices  = [
        1 - melSimOrthoDf["mel_min_sig_ratio_avg_f_m_apn_uq_ff"],
        melSimOrthoDf["mel_min_sig_ratio_avg_m_f_apn_uq_ff"] - 1,
        1 - melSimOrthoDf["mel_mean_ratio_avg_f_m_apn_uq_ff"],
        melSimOrthoDf["mel_mean_ratio_avg_m_f_apn_uq_ff"] - 1,
        1,
        -1
    ]
    melSimOrthoDf["mel_ratio"] = np.select(melRatioConditions, melRatioChoices, "oops").astype(float)

    # Sim ratios
    simRatioConditions = [
        melSimOrthoDf["flag_M_sim"] == 1,
        melSimOrthoDf["flag_F_sim"] == 1,
        (melSimOrthoDf["flag_U_sim"] == 1)
            & (melSimOrthoDf["sim_mean_ratio_avg_f_m_apn_uq_ff"] < 1)        
            & (melSimOrthoDf["sim_mean_ratio_avg_m_f_apn_uq_ff"] > 1),
        (melSimOrthoDf["flag_U_sim"] == 1)
            & (melSimOrthoDf["sim_mean_ratio_avg_m_f_apn_uq_ff"] < 1)
            & (melSimOrthoDf["sim_mean_ratio_avg_f_m_apn_uq_ff"] > 1),
        (melSimOrthoDf["flag_U_sim"] == 1)
            & (melSimOrthoDf["sim_mean_ratio_avg_m_f_apn_uq_ff"] > melSimOrthoDf["sim_mean_ratio_avg_f_m_apn_uq_ff"]),
        (melSimOrthoDf["flag_U_sim"] == 1)
            & (melSimOrthoDf["sim_mean_ratio_avg_m_f_apn_uq_ff"] < melSimOrthoDf["sim_mean_ratio_avg_f_m_apn_uq_ff"]),
    ]
    simRatioChoices  = [
        1 - melSimOrthoDf["sim_min_sig_ratio_avg_f_m_apn_uq_ff"],
        melSimOrthoDf["sim_min_sig_ratio_avg_m_f_apn_uq_ff"] - 1,
        1 - melSimOrthoDf["sim_mean_ratio_avg_f_m_apn_uq_ff"],
        melSimOrthoDf["sim_mean_ratio_avg_m_f_apn_uq_ff"] - 1,
        1,
        -1
    ]
    melSimOrthoDf["sim_ratio"] = np.select(simRatioConditions, simRatioChoices, "oops").astype(float)

    # Output divergent sex biased genes
    melSimOrthoDf[[
        "mel_gene_id",
        "mel_geneSymbol",
        "sim_gene_id",
        "sim_geneSymbol",
        "mel_ratio",
        "sim_ratio"] + GOflagLst + expFlagLst + mappingFlagLst].to_csv(args.outDir + "/ortho_MF_ratio_all_gene.csv", index=False)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
