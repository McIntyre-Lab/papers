#!/usr/bin/env python

import sss_upset_plot_functions_01ksb as plotter
import glob
import pandas as pd


def main():

    # Define ind and outd
    indir = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/"
    outdir = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/manuscript"

    # Create jxnHash upset plot for all 5 genomes
    mergeFlagFiles = glob.glob(f"{indir}/flag_fiveSpecies_2_*_ujc.csv")
    for flagFile in mergeFlagFiles:
        plotter.plotSharedJxnHashPerGenome(flagFile=flagFile, outdir=outdir)

    # Create jxnHash upset plot for below 4 genes on *ALL GENOME*
    # 5 Gene == Fru (FBgn0004652), Dsx (FBgn0000504), Sxl (FBgn0264270), Doa (FBgn0264270), Tra2 (FBgn0003742)
    fiveGeneDct = {
        "Doa": "FBgn0265998",
        "Dsx": "FBgn0000504",
        "Fru": "FBgn0004652",
        "Sxl": "FBgn0264270",
        "Tra2": "FBgn0003742"
    }

    for genome, otherName in {"dmel6": None, "dsim2": 'dsim', "dsan1": 'dsan', "dyak2": 'dyak', "dser1": 'dser'}.items():

        if genome == 'dmel6':
            geneDct = fiveGeneDct
        else:
            speciesLnkFile = f"{indir}/{otherName}_sexdet.csv"
            speciesLnkDf = pd.read_csv(speciesLnkFile, low_memory=False)
            speciesLnkDf = speciesLnkDf[[
                col for col in speciesLnkDf.columns if 'geneID' in col]]

        plotter.plotGeneJxnHash(
            f"{indir}/fiveSpecies_2_{genome}_ujc_gene_key.csv",
            f"{indir}/flag_fiveSpecies_2_{genome}_ujc.csv",
            geneDct,
            genome,
            outdir
        )

    for genome in ["dmel6", "dsim2", "dsan1", "dyak2", "dser1"]:
        flagFile = ("/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/"
                    "sex_specific_splicing/cross_species_link_files/"
                    f"flag_fiveSpecies_2_{genome}_ujc.csv")

        erpFile = ("/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/"
                   "sex_specific_splicing/cross_species_link_files/"
                   f"fiveSpecies_2_{genome}_ujc_er_vs_fiveSpecies_2_{genome}_ujc_ERP.csv")

        plotter.plotSharedERPPerGenome(
            flagFile,
            erpFile,
            genome,
            outdir
        )


if __name__ == '__main__':
    main()
