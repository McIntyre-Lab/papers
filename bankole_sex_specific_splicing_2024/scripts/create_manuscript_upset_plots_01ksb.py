#!/usr/bin/env python

import sss_upset_plot_functions_01ksb as plotter
import glob


def main():

    # Define ind and outd
    indir = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/cross_species_link_files/"
    outdir = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/manuscript"

    # Create jxnHash upset plot for all 5 genomes
    mergeFlagFiles = glob.glob(f"{indir}/flag_fiveSpecies_2_*_ujc.csv")
    for flagFile in mergeFlagFiles:
        plotter.plotSharedJxnHashPerGenome(flagFile=flagFile, outdir=outdir)

    # Create jxnHash upset plot for below 4 genes on *MEL*
    # 4 Gene == Fru (FBgn0004652), Dsx (FBgn0000504), Sxl (FBgn0264270), Doa (FBgn0264270)
    fourGeneDct = {
        "Doa": "FBgn0265998",
        "Dsx": "FBgn0000504",
        "Fru": "FBgn0004652",
        "Sxl": "FBgn0264270"
    }
    genome = "dmel6"

    plotter.plotGeneJxnHash(
        f"{indir}/fiveSpecies_2_dmel6_ujc_gene_key.csv",
        f"{indir}/flag_fiveSpecies_2_dmel6_ujc.csv",
        fourGeneDct,
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
