#!/usr/bin/env python

import sss_upset_plot_functions_01ksb as plotter
import glob
import pandas as pd


def main():

    # Define ind and outd
    indir = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/"
    outdir = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Figures/Upset"

    # Create jxnHash upset plot for all 5 genomes
    flagFileLst = glob.glob(f"{indir}/flag_fiveSpecies_2_*_ujc.csv")
    for flagFile in flagFileLst:
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

        geneDct = dict()

        if genome == 'dmel6':
            geneDct = fiveGeneDct
        else:
            speciesLnkFile = f"{indir}/{otherName}_sexdet_w_gtf_info.csv"
            speciesLnkDf = pd.read_csv(speciesLnkFile, low_memory=False)

            for row in speciesLnkDf.to_dict('records'):
                if row[f'num_{genome}_geneID'] == 1:
                    if row['mel_symbol'] == "Doa":
                        geneDct['Doa'] = row[f'{genome}_geneID']
                    elif row['mel_symbol'] == "dsx":
                        geneDct['Dsx'] = row[f'{genome}_geneID']
                    elif row['mel_symbol'] == "fru":
                        geneDct['Fru'] = row[f'{genome}_geneID']
                    elif row['mel_symbol'] == "Sxl":
                        geneDct['Sxl'] = row[f'{genome}_geneID']
                    if row['mel_symbol'] == "tra2":
                        geneDct['Tra2'] = row[f'{genome}_geneID']

            # Manually set the two non-1to1 geneIDs
            if genome == 'dsim2':
                geneDct['Fru'] = "FBgn0043395"
            if genome == 'dser1':
                geneDct['Fru'] = "LOC110187623"
                geneDct['Doa'] = "LOC110184797"

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
