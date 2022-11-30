#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=(
        "Prepare datasets for scatter plots of gene-level expression values from "
        "alignments of D. melanogaster vs. D. simulans to the same genome."
        )
    )

    # Input data
    parser.add_argument(
        "-m",
        "--mel",
        dest="mel",
        required=True,
        help="CSV of gene-level expression values for D. melanogaster and D. simulans on D. melanogaster coordinates."
    )
    parser.add_argument(
        "-sf",
        "--sim",
        dest="sim",
        required=True,
        help="CSV of gene-level expression values for D. melanogaster and D. simulans on D. simulans coordinates."
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
    melCoordDf = pd.read_csv(args.mel, low_memory=False)
    simCoordDf = pd.read_csv(args.sim, low_memory=False)
    orthoDf = pd.read_csv(args.ortho, low_memory=False)

    # melCoordDf = pd.read_csv("/Users/adalena/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq/Coverage_counts/gene_cvrg_cnts_ortho/cvrGene_ortho2mel_log_uq_apn.csv", low_memory=False)
    # simCoordDf = pd.read_csv("/Users/adalena/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq/Coverage_counts/gene_cvrg_cnts_ortho/cvrGene_ortho2sim_log_uq_apn.csv", low_memory=False)
    # orthoDf = pd.read_csv("/Users/adalena/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/supp_files/dmel_dsim_ortholog_chip_rna_flags.csv", low_memory=False)

    # Get flags from ortholog gene-level results file
    GOflagLst = ['biological_process', 'goterm_biol_process', 'molecular_function',
                 'goterm_mol_function', 'cellular_component', 'goterm_cell_component',]
    chromFlagLst = ['flag_has_male_k4_mel', 'flag_has_male_k4_sim', 'flag_has_female_k4_mel',
                    'flag_has_female_k4_sim', 'flag_has_male_k27_mel', 'flag_has_male_k27_sim',
                    'flag_has_female_k27_mel', 'flag_has_female_k27_sim']

    # Split mel and sim gene_id values in gene expression files
    melCoordDf["mel_geneID"] = melCoordDf["melID_simID"].str.split("_").str[0]
    melCoordDf["sim_geneID"] = melCoordDf["melID_simID"].str.split("_").str[1]
    simCoordDf["mel_geneID"] = simCoordDf["melID_simID"].str.split("_").str[0]
    simCoordDf["sim_geneID"] = simCoordDf["melID_simID"].str.split("_").str[1]

    # Get average male and female expression for each species on each set of coordinates
    melCoordDf["avg_m_mel"] = melCoordDf[[c for c in melCoordDf.columns if "mel_" in c and "_m_" in c]].mean(axis=1)
    melCoordDf["avg_f_mel"] = melCoordDf[[c for c in melCoordDf.columns if "mel_" in c and "_f_" in c]].mean(axis=1)
    melCoordDf["avg_m_sim"] = melCoordDf[[c for c in melCoordDf.columns if "sim_" in c and "_m_" in c]].mean(axis=1)
    melCoordDf["avg_f_sim"] = melCoordDf[[c for c in melCoordDf.columns if "sim_" in c and "_f_" in c]].mean(axis=1)
    simCoordDf["avg_m_mel"] = simCoordDf[[c for c in simCoordDf.columns if "mel_" in c and "_m_" in c]].mean(axis=1)
    simCoordDf["avg_f_mel"] = simCoordDf[[c for c in simCoordDf.columns if "mel_" in c and "_f_" in c]].mean(axis=1)
    simCoordDf["avg_m_sim"] = simCoordDf[[c for c in simCoordDf.columns if "sim_" in c and "_m_" in c]].mean(axis=1)
    simCoordDf["avg_f_sim"] = simCoordDf[[c for c in simCoordDf.columns if "sim_" in c and "_f_" in c]].mean(axis=1)

    # Merge each set of coordinates with the ortholog results by mel_geneID and sim_geneID to get gene groups
    melMergeDf = pd.merge(
        melCoordDf[["mel_geneID", "sim_geneID", "avg_m_mel", "avg_f_mel", "avg_m_sim", "avg_f_sim"]],
        orthoDf,
        on=["mel_geneID", "sim_geneID"],
        how="outer",
        validate="1:1",
        indicator="merge_check"
    )
    # NOTE: There are 4 genes in the gene expression table that are NOT in the ortholog results
    # What are these?
    # FBgn0024232  FBgn0269679 - the sim gene is missing from the table because it has all multigene fragments
    melMergeDf = melMergeDf[melMergeDf["merge_check"]=="both"].drop(columns=["merge_check"])

    simMergeDf = pd.merge(
        simCoordDf[["mel_geneID", "sim_geneID", "avg_m_mel", "avg_f_mel", "avg_m_sim", "avg_f_sim"]],
        orthoDf,
        on=["mel_geneID", "sim_geneID"],
        how="outer",
        validate="1:1",
        indicator="merge_check"
    )
    # NOTE: There are 123 genes in the gene expression table that are NOT in the ortholog results
    # What are these?
    simMergeDf = simMergeDf[simMergeDf["merge_check"]=="both"].drop(columns=["merge_check"])

    # Merge the coordinates and output the mel on mel vs. sim and the sim on mel vs. sim
    bothCoordDf = pd.merge(
        melMergeDf,
        simCoordDf[["mel_geneID", "sim_geneID", "avg_m_mel", "avg_f_mel", "avg_m_sim", "avg_f_sim"]],
        on=["mel_geneID", "sim_geneID"],
        suffixes=["2mel", "2sim"],
        how="outer",
        validate="1:1",
        indicator="merge_check"
    )
    # bothCoordDf["merge_check"].value_counts()
        # both          9328
        # right_only     478
        # left_only      145
    bothCoordDf = bothCoordDf[bothCoordDf["merge_check"]=="both"].drop(columns=["merge_check"])

    # Make ratios
    bothCoordDf["ratio_mel2mel"] = np.where(
        bothCoordDf["avg_f_mel2mel"].abs() > bothCoordDf["avg_m_mel2mel"].abs(),
        1 - (bothCoordDf["avg_m_mel2mel"]/bothCoordDf["avg_f_mel2mel"]).abs(),
        1 - (bothCoordDf["avg_f_mel2mel"]/bothCoordDf["avg_m_mel2mel"]).abs()
    )
    bothCoordDf["ratio_sim2mel"] = np.where(
        bothCoordDf["avg_f_sim2mel"].abs() > bothCoordDf["avg_m_sim2mel"].abs(),
        1 - (bothCoordDf["avg_m_sim2mel"]/bothCoordDf["avg_f_sim2mel"]).abs(),
        1 - (bothCoordDf["avg_f_sim2mel"]/bothCoordDf["avg_m_sim2mel"]).abs()
    )
    bothCoordDf["ratio_mel2sim"] = np.where(
        bothCoordDf["avg_f_mel2sim"].abs() > bothCoordDf["avg_m_mel2sim"].abs(),
        1 - (bothCoordDf["avg_m_mel2sim"]/bothCoordDf["avg_f_mel2sim"]).abs(),
        1 - (bothCoordDf["avg_f_mel2sim"]/bothCoordDf["avg_m_mel2sim"]).abs()
    )
    bothCoordDf["ratio_sim2sim"] = np.where(
        bothCoordDf["avg_f_sim2sim"].abs() > bothCoordDf["avg_m_sim2sim"].abs(),
        1 - (bothCoordDf["avg_m_sim2sim"]/bothCoordDf["avg_f_sim2sim"]).abs(),
        1 - (bothCoordDf["avg_f_sim2sim"]/bothCoordDf["avg_m_sim2sim"]).abs()
    )

    # Flag genes that map at least 2x better in one set of coord vs the other in males and females
    bothCoordDf["flag_map_better_f_mel2mel"] = np.where(
        np.exp(bothCoordDf["avg_f_mel2mel"]) >= 2 * np.exp(bothCoordDf["avg_f_mel2sim"]),
        1,
        0
    )
    bothCoordDf["flag_map_better_m_mel2mel"] = np.where(
        np.exp(bothCoordDf["avg_m_mel2mel"]) >= 2 * np.exp(bothCoordDf["avg_m_mel2sim"]),
        1,
        0
    )
    bothCoordDf["flag_map_better_f_sim2sim"] = np.where(
        np.exp(bothCoordDf["avg_f_sim2sim"]) >= 2 * np.exp(bothCoordDf["avg_f_sim2mel"]),
        1,
        0
    )
    bothCoordDf["flag_map_better_m_sim2sim"] = np.where(
        np.exp(bothCoordDf["avg_m_sim2sim"]) >= 2 * np.exp(bothCoordDf["avg_m_sim2mel"]),
        1,
        0
    )
    bothCoordDf["flag_map_better_f_mel2sim"] = np.where(
        np.exp(bothCoordDf["avg_f_mel2sim"]) >= 2 * np.exp(bothCoordDf["avg_f_mel2mel"]),
        1,
        0
    )
    bothCoordDf["flag_map_better_m_mel2sim"] = np.where(
        np.exp(bothCoordDf["avg_m_mel2sim"]) >= 2 * np.exp(bothCoordDf["avg_m_mel2mel"]),
        1,
        0
    )
    bothCoordDf["flag_map_better_f_sim2mel"] = np.where(
        np.exp(bothCoordDf["avg_f_sim2mel"]) >= 2 * np.exp(bothCoordDf["avg_f_sim2sim"]),
        1,
        0
    )
    bothCoordDf["flag_map_better_m_sim2mel"] = np.where(
        np.exp(bothCoordDf["avg_m_sim2mel"]) >= 2 * np.exp(bothCoordDf["avg_m_sim2sim"]),
        1,
        0
    )

    # Loop over coordinate sets an output conserved and diverged datasets
    for coord in ["mel", "sim", "bothCoord"]:
        if coord == "mel":
            coordMergeDf = melMergeDf.copy()
            expFlagLst = ["avg_m_mel", "avg_f_mel", "avg_m_sim", "avg_f_sim"]
        elif coord == "sim":
            coordMergeDf = simMergeDf.copy()
            expFlagLst = ["avg_m_mel", "avg_f_mel", "avg_m_sim", "avg_f_sim"]
        else:
            coordMergeDf = bothCoordDf.copy()
            expFlagLst = ["avg_m_mel2mel", "avg_f_mel2mel", "avg_m_sim2mel", "avg_f_sim2mel",
                          "avg_m_mel2sim", "avg_f_mel2sim", "avg_m_sim2sim", "avg_f_sim2sim",
                          "ratio_mel2mel", "ratio_sim2mel",
                          "ratio_mel2sim", "ratio_sim2sim",
                          "flag_map_better_f_mel2mel", "flag_map_better_m_mel2mel",
                          "flag_map_better_f_sim2sim", "flag_map_better_m_sim2sim",
                          "flag_map_better_f_mel2sim", "flag_map_better_m_mel2sim",
                          "flag_map_better_f_sim2mel", "flag_map_better_m_sim2mel",
                          "flag_conserved_male_bias", "flag_conserved_female_bias"]
            # Export all for both coords
            coordMergeDf[[
                "mel_geneID",
                "mel_geneSymbol",
                "sim_geneID",
                "sim_geneSymbol"] + expFlagLst + GOflagLst + chromFlagLst].to_csv(
                    args.outDir + "/ortho_map2"+coord+"_all.csv", index=False)

        # Get conserved male-biased and female-biased genes
        MconservDf = coordMergeDf[
            (coordMergeDf["flag_conserved_male_bias"]==1)
            & (coordMergeDf["xsome"].isin(["X", "A"]))
        ].copy()
        FconservDf = coordMergeDf[
            (coordMergeDf["flag_conserved_female_bias"]==1)
            & (coordMergeDf["xsome"].isin(["X", "A"]))
        ].copy()
    
        # Output convserved ratio CSV
        MconservDf[[
            "mel_geneID",
            "mel_geneSymbol",
            "sim_geneID",
            "sim_geneSymbol"] + expFlagLst + GOflagLst + chromFlagLst].to_csv(
                args.outDir + "/ortho_map2"+coord+"_conserved_M.csv", index=False)
        FconservDf[[
            "mel_geneID",
            "mel_geneSymbol",
            "sim_geneID",
            "sim_geneSymbol"] + expFlagLst + GOflagLst + chromFlagLst].to_csv(
                args.outDir + "/ortho_map2"+coord+"_conserved_F.csv", index=False)


        # Get divergent sex biased genes
        # MmelFsimDf = coordMergeDf[
        #     (coordMergeDf["flag_M_mel_F_sim"]==1)
        #     & (coordMergeDf["xsome"].isin(["X", "A"]))
        # ].copy()

        # FmelMsimDf = coordMergeDf[
        #     (coordMergeDf["flag_M_sim_F_mel"]==1)
        #     & (coordMergeDf["xsome"].isin(["X", "A"]))
        # ].copy()

        # # Output divergent sex biased genes
        # MmelFsimDf[[
        #     "mel_geneID",
        #     "mel_geneSymbol",
        #     "sim_geneID",
        #     "sim_geneSymbol",
        #     "avg_m_mel",
        #     "avg_f_mel",
        #     "avg_m_sim",
        #     "avg_f_sim"] + GOflagLst + chromFlagLst].to_csv(args.outDir + "/ortho_map2"+coord+"_M_mel_F_sim.csv", index=False)
        # FmelMsimDf[[
        #     "mel_geneID",
        #     "mel_geneSymbol",
        #     "sim_geneID",
        #     "sim_geneSymbol",
        #     "avg_m_mel",
        #     "avg_f_mel",
        #     "avg_m_sim",
        #     "avg_f_sim"] + GOflagLst + chromFlagLst].to_csv(args.outDir + "/ortho_map2"+coord+"_F_mel_M_sim.csv", index=False)

    


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
