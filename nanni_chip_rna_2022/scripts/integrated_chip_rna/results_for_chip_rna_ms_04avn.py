#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import mcnemar
from scipy.stats import binom_test

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Run tests for results in Drosophila ChIP-RNA manuscript."
    )

    # Input data
    parser.add_argument(
        "-mg",
        "--mel-gene",
        dest="inMG",
        required=True,
        help=(
            "Input CSV of D. melanogaster gene-level file (from supplement)."
        )
    )
    parser.add_argument(
        "-mf",
        "--mel-feature",
        dest="inMF",
        required=True,
        help=(
            "Input CSV of D. melanogaster feature-level file (from supplement)."
        )
    )
    parser.add_argument(
        "-sg",
        "--sim-gene",
        dest="inSG",
        required=True,
        help=(
            "Input CSV of D. simulans gene-level file (from supplement)."
        )
    )
    parser.add_argument(
        "-sf",
        "--sim-feature",
        dest="inSF",
        required=True,
        help=(
            "Input CSV of D. simulans feature-level file (from supplement)."
        )
    )
    parser.add_argument(
        "-or",
        "--ortho",
        dest="inO",
        required=True,
        help=(
            "Input CSV of D. melangoster and D. simulans orthologous "
            "gene-level file (from supplement)."
        )
    )
    parser.add_argument(
        "-b",
        "--map-bias",
        dest="inMapBias",
        required=False,
        action="store_true",
        help=(
            "Exclude all genes with any mapping bias in males or females of "
            "either species (exclude where flag_map_better_f_mel2mel + flag_map_better_m_mel2mel + "
            "flag_map_better_f_sim2sim + flag_map_better_m_sim2sim + flag_map_better_f_mel2sim "
            "+ flag_map_better_m_mel2sim + flag_map_better_f_sim2mel + flag_map_better_m_sim2mel > 0)."
            " NOTE: Feature tests are not performed."
        )
    )
    parser.add_argument(
        "-c",
        "--consistent-map-bias",
        dest="inConsMapBias",
        required=False,
        action="store_true",
        help=(
            "Exclude genes with a consistent mapping bias in males and females of "
            "both species (exclude where flag_map_better_2_mel + flag_map_better_2_sim > 0)."
            " NOTE: Feature tests are not performed."
        )
    )

    # Output data
    parser.add_argument(
        "-d",
        "--output-dir",
        dest="outDir",
        required=True,
        help="Output directory for output files of test results."
    )

    args = parser.parse_args()
    return args

def indent_df(data, num_indent=1):
    spaces = " " * (num_indent*8)
    return spaces + data.to_string().replace("\n", "\n"+spaces)

def test_XA(data, flag, file, title):
    # Get crosstab of given flag vs X/A within the given dataframe
    crosstabDF = pd.crosstab(data[flag], data["xsome"])[["X", "A"]].T
    # Calculate chi2
    chiSquare = chi2_contingency(crosstabDF)
    # Calculate Cramer's V
    cramersV = np.sqrt((chiSquare[0]/crosstabDF.sum().sum()))
    # Calculate Fisher exact and odds ratio
    fisherExact2 = fisher_exact(crosstabDF, alternative="two-sided")
    # test that the first value in the table is less than by chance (H0: >=)
    fisherExactLess = fisher_exact(crosstabDF, alternative="less")
    # test that the first value in the table is greater than by chance (H0: <=)
    fisherExactGreater = fisher_exact(crosstabDF, alternative="greater")
    # Calculate McNemar test of homogeneity
    MN = mcnemar(crosstabDF.to_numpy())
#    return crosstabDF, chiSquare, cramersV, fisherExact2, fisherExactLT, fisherExactGT
    print_test_results(file, title, crosstabDF, chiSquare, cramersV, fisherExact2, fisherExactLess, fisherExactGreater, MN)

def test_2var(data, flag1, flag2, file, title):
    # Get crosstab of given flag vs X/A within the given dataframe
    crosstabDF = pd.crosstab(data[flag1], data[flag2])
    # Calculate chi2 test of independence
    chiSquare = chi2_contingency(crosstabDF)
    # Calculate Cramer's V
    cramersV = np.sqrt((chiSquare[0]/crosstabDF.sum().sum()))
    # Calculate Fisher exact and odds ratio
    fisherExact2 = fisher_exact(crosstabDF, alternative='two-sided')
    # test that the first value in the table is less than by chance (H0: >=)
    fisherExactLess = fisher_exact(crosstabDF, alternative='less')
    # test that the first value in the table is greater than by chance (H0: <=)
    fisherExactGreater = fisher_exact(crosstabDF, alternative='greater')
    # Binomial test two tailed
    binomial2 = binom_test(
        x=max(crosstabDF[1][0], crosstabDF[0][1]),
        n=crosstabDF[1][0]+crosstabDF[0][1],
        p=0.5,
        alternative='two-sided'
    )
    # Calculate McNemar test of homogeneity
    MN = mcnemar(crosstabDF.to_numpy())
    # Calculate Kappa = (po - pe) / (1 - pe)
    #   po = (a + d) / total = proportion of diagonal where flags are the same
    #   pe = pyes + pno
    #       pyes = ((a+b) / total) * ((a+c) / total)
    #           = prop of first row * prop of first column
    #       pno = ((c+d) / total) * ((b+d) / total)
    #           = prop of second row * prop of second column
    total = crosstabDF.sum().sum()
    po = (crosstabDF[0][0] + crosstabDF[1][1]) / total
    pyes = (
        (
            (crosstabDF[0][0] + crosstabDF[1][0]) / total
        ) * (
            (crosstabDF[0][0] + crosstabDF[0][1]) / total
        )
    )
    pno = (
        (
            (crosstabDF[0][1] + crosstabDF[1][1]) / total
        ) * (
            (crosstabDF[1][0] + crosstabDF[1][1]) / total
        )
    )
    pe = pyes + pno
    kappa = (po - pe) / (1 - pe)
    print_test_results(file, title, crosstabDF, chiSquare, cramersV, fisherExact2, fisherExactLess, fisherExactGreater, MN, kappa, binomial2)

def print_test_results(file, title, crosstabDF, chiSquare, cramersV, fisherExact2, fisherExactLess, fisherExactGreater, MN, kappa = None, binomial2 = None):
    outText = (
        "*** {} ***\n\tObserved Counts:\n{}\n\n\tExpected Counts:\n{}\n\n"
        "\tChi2 = {}, pval = {}, Cramer's V = {}\n"
        "\tOdds Ratio = {}\n"
        "\tFisher exact pval (two-sided) = {}\n"
        "\tFisher exact pval (less) = {}\n"
        "\tFisher exact pval (greater) = {}\n"        
        "\tMcNemar = {}, pval = {}\n".format(
                title,
                indent_df(crosstabDF),
                indent_df(pd.DataFrame(
                        chiSquare[3],
                        index=crosstabDF.index,
                        columns=crosstabDF.columns
                    )),
                chiSquare[0],
                chiSquare[1],
                cramersV,
                fisherExact2[0],
                fisherExact2[1],
                fisherExactLess[1],
                fisherExactGreater[1],
                MN.statistic,
                MN.pvalue
        )
    )
    if kappa is not None:
        outText = outText + "\tKappa = {}\n".format(kappa)
    if binomial2 is not None:
        outText = outText + "\tP(X = {}), ~Binomial(n = {}+{} = {}, p = 0.5) = {}\n\n\n".format(
            max(crosstabDF[1][0], crosstabDF[0][1]),
            crosstabDF[1][0],
            crosstabDF[0][1],
            crosstabDF[1][0] + crosstabDF[0][1],
            binomial2
        )
    else:
        outText = outText + "\n\n"
    file.write(outText)


def main():
    # Get gene inputs
    melDF = pd.read_csv(args.inMG, low_memory=False)
    simDF = pd.read_csv(args.inSG, low_memory=False)
    orthoAllDF = pd.read_csv(args.inO, low_memory=False)

    if args.inMapBias and args.inConsMapBias:
        print("!!!ERROR: Must provide eith -b or -c, not both.")
        exit()
    if args.inMapBias:
        melDF = melDF[melDF[[c for c in melDF.columns if "flag_map_better_" in c and "_2_" not in c]].fillna(0).sum(axis=1)==0].copy()
        simDF = simDF[simDF[[c for c in simDF.columns if "flag_map_better_" in c and "_2_" not in c]].fillna(0).sum(axis=1)==0].copy()
        orthoAllDF = orthoAllDF[orthoAllDF[[c for c in orthoAllDF.columns if "flag_map_better_" in c and "_2_" not in c]].fillna(0).sum(axis=1)==0].copy()
    if args.inConsMapBias:
        melDF = melDF[melDF[[c for c in melDF.columns if "flag_map_better_2_" in c]].sum(axis=1)==0].copy()
        simDF = simDF[simDF[[c for c in simDF.columns if "flag_map_better_2_" in c]].sum(axis=1)==0].copy()
        orthoAllDF = orthoAllDF[orthoAllDF[[c for c in orthoAllDF.columns if "flag_map_better_2_" in c]].sum(axis=1)==0].copy()

    if not args.inMapBias and not args.inConsMapBias:
        # Get feature inputs
        melFeature = pd.read_csv(args.inMF, low_memory=False)
        simFeature = pd.read_csv(args.inSF, low_memory=False)
        # Add flags to feature input - !!! Move this to different script for feature flags
        melFeature["flag_any_K4_on"] = np.where(
                (melFeature["flag_f_K4_on"] == 1)
                | (melFeature["flag_m_K4_on"] == 1),
                1,
                0
        )
        simFeature["flag_any_K4_on"] = np.where(
                (simFeature["flag_f_K4_on"] == 1)
                | (simFeature["flag_m_K4_on"] == 1),
                1,
                0
        )
        melFeature["flag_any_K27_on"] = np.where(
                (melFeature["flag_f_K27_on"] == 1)
                | (melFeature["flag_m_K27_on"] == 1),
                1,
                0
        )
        simFeature["flag_any_K27_on"] = np.where(
                (simFeature["flag_f_K27_on"] == 1)
                | (simFeature["flag_m_K27_on"] == 1),
                1,
                0
        )

    # Drop orthologs pairs that are not one-2-one
    orthoCount = open("{}/ortholog_counts.txt".format(args.outDir),"w")
    orthoCount.write("\n{} total orthologous gene pairs between D. melanogsater and D. simulans\n".format(
        len(orthoAllDF)
    ))

    orthoDF = orthoAllDF[orthoAllDF["flag_one2one_ortholog"]==1].copy()
    orthoDF["mel_xsome_scaffold"] = orthoDF["xsome"].fillna("Scaffold")
    orthoDF["sim_xsome_scaffold"] = orthoDF["sim_xsome"].fillna("Scaffold")

    orthoCount.write("\t{} one-to-one ortholgous gene pairs\n".format(len(orthoDF)))
    orthoCount.write("\nChromosome pairs of one-to-one ortholgous genes:\n{}\n".format(
        pd.crosstab(orthoDF["mel_xsome_scaffold"], orthoDF["sim_xsome_scaffold"])))

    # Drop ortholog pairs that the chromosome does not match (e.g. X in mel but A in sim)
    orthoDF = orthoDF[orthoDF["xsome"] == orthoDF["sim_xsome"]]

    # Get orthologs pairs that the chromosomse are the same and on X or A
    orthoXA = orthoDF[orthoDF["xsome"].isin(["X", "A"])].copy()
    orthoCount.write(
        "\n{} one-to-one orthologous genes are on the X and autosomes of both "
        "species (both on X or both on A)\n".format(
        len(orthoXA)
    ))
    orthoCount.close()
    melOrthoGenes = orthoXA["mel_geneID"].unique()
    simOrthoGenes = orthoXA["sim_geneID"].unique()

    # Drop genes not expressed
    melExpressWsexLim = melDF[melDF["flag_expressed"] == 1].copy()
    simExpressWsexLim = simDF[simDF["flag_expressed"] == 1].copy()
    orthoExpressBothWsexLim = orthoDF[
            (orthoDF["flag_expressed_mel"]==1)
            & (orthoDF["flag_expressed_sim"]==1)
    ].copy()
    orthoExpressEitherWsexLim = orthoDF[
            (orthoDF["flag_expressed_mel"]==1)
            | (orthoDF["flag_expressed_sim"]==1)
    ].copy()
    orthoExpressEitherXA = orthoExpressEitherWsexLim[orthoExpressEitherWsexLim["xsome"].isin(["X", "A"])]
    
    # Drop sex-limited genes
    melExpress = melExpressWsexLim[melExpressWsexLim["flag_sex_limited"] == 0]
    simExpress = simExpressWsexLim[simExpressWsexLim["flag_sex_limited"] == 0]
    orthoExpress = orthoExpressBothWsexLim[
            (orthoExpressBothWsexLim["flag_sex_limited_mel"]==0)
            & (orthoExpressBothWsexLim["flag_sex_limited_sim"]==0)
    ]

    # Get dataframe subsets of X and A (no chrom 4) only
    melXA = melDF[melDF["xsome"].isin(["X", "A"])]
    simXA = simDF[simDF["xsome"].isin(["X", "A"])]
    orthoXA = orthoDF[orthoDF["xsome"].isin(["X", "A"])]
    melExpressXA = melExpress[melExpress["xsome"].isin(["X", "A"])]
    simExpressXA = simExpress[simExpress["xsome"].isin(["X", "A"])]
    orthoExpressXA = orthoExpress[orthoExpress["xsome"].isin(["X", "A"])].copy()

    # Get counts for Table 1
    table1Count = open("{}/Table_1_counts.tsv".format(args.outDir),"w")
    table1Count.write("\tD. melanogaster (X,A)\tD. simulans (X,A)\tOrthologs (X,A)\t% Orthologs Mel Male Open\t% Orthologs Mel Female Open\t% Orthologs Mel Male Closed\t% Orthologs Mel Female Closed\t% Orthologs Sim Male Open\t% Orthologs Sim Female Open\t% Orthologs Sim Male Closed\t% Orthologs Sim Female Closed\n")
    table1Count.write("Male-biased\t{0} ({1}, {2})\t{3} ({4}, {5})\t{6} ({7}, {8})\t{9:.2%}\t{10:.2%}\t{11:.2%}\t{12:.2%}\t{13:.2%}\t{14:.2%}\t{15:.2%}\t{16:.2%}\n".format(
        melExpressXA["flag_M"].sum(),
        melExpressXA[melExpressXA["xsome"]=="X"]["flag_M"].sum(),
        melExpressXA[melExpressXA["xsome"]=="A"]["flag_M"].sum(),
        simExpressXA["flag_M"].sum(),
        simExpressXA[simExpressXA["xsome"]=="X"]["flag_M"].sum(),
        simExpressXA[simExpressXA["xsome"]=="A"]["flag_M"].sum(),
        orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_conserved_male_bias"].sum() / orthoExpressXA["flag_conserved_male_bias"].sum(),
    ))
    table1Count.write("Female-biased\t{0} ({1}, {2})\t{3} ({4}, {5})\t{6} ({7}, {8})\t{9:.2%}\t{10:.2%}\t{11:.2%}\t{12:.2%}\t{13:.2%}\t{14:.2%}\t{15:.2%}\t{16:.2%}\n".format(
        melExpressXA["flag_F"].sum(),
        melExpressXA[melExpressXA["xsome"]=="X"]["flag_F"].sum(),
        melExpressXA[melExpressXA["xsome"]=="A"]["flag_F"].sum(),
        simExpressXA["flag_F"].sum(),
        simExpressXA[simExpressXA["xsome"]=="X"]["flag_F"].sum(),
        simExpressXA[simExpressXA["xsome"]=="A"]["flag_F"].sum(),
        orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_conserved_female_bias"].sum() / orthoExpressXA["flag_conserved_female_bias"].sum(),
    ))
    table1Count.write("Male- and Female-biased\t{0} ({1}, {2})\t{3} ({4}, {5})\t{6} ({7}, {8})\t{9:.2%}\t{10:.2%}\t{11:.2%}\t{12:.2%}\t{13:.2%}\t{14:.2%}\t{15:.2%}\t{16:.2%}\n".format(
        melExpressXA["flag_MF"].sum(),
        melExpressXA[melExpressXA["xsome"]=="X"]["flag_MF"].sum(),
        melExpressXA[melExpressXA["xsome"]=="A"]["flag_MF"].sum(),
        simExpressXA["flag_MF"].sum(),
        simExpressXA[simExpressXA["xsome"]=="X"]["flag_MF"].sum(),
        simExpressXA[simExpressXA["xsome"]=="A"]["flag_MF"].sum(),
        orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_conserved_male_and_female"].sum() / orthoExpressXA["flag_conserved_male_and_female"].sum(),
    ))
    table1Count.write("Unbiased\t{0} ({1}, {2})\t{3} ({4}, {5})\t{6} ({7}, {8})\t{9:.2%}\t{10:.2%}\t{11:.2%}\t{12:.2%}\t{13:.2%}\t{14:.2%}\t{15:.2%}\t{16:.2%}\n".format(
        melExpressXA["flag_U"].sum(),
        melExpressXA[melExpressXA["xsome"]=="X"]["flag_U"].sum(),
        melExpressXA[melExpressXA["xsome"]=="A"]["flag_U"].sum(),
        simExpressXA["flag_U"].sum(),
        simExpressXA[simExpressXA["xsome"]=="X"]["flag_U"].sum(),
        simExpressXA[simExpressXA["xsome"]=="A"]["flag_U"].sum(),
        orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_conserved_unbiased"].sum() / orthoExpressXA["flag_conserved_unbiased"].sum(),
    ))
    table1Count.write("Reversal\tMale\tFemale\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_M_mel_F_sim"].sum() / orthoExpressXA["flag_M_mel_F_sim"].sum(),
    ))
    table1Count.write("Reversal\tFemale\tMale\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_M_sim_F_mel"].sum() / orthoExpressXA["flag_M_sim_F_mel"].sum(),
    ))
    table1Count.write("Gain/Loss\tMale\tMale and Female\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)&(orthoExpressXA["xsome"]=="X")]),
        len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)&(orthoExpressXA["xsome"]=="A")]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_mel"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_mel"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_mel"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_mel"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_sim"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_sim"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_sim"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_sim"]==1)&(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_M_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
    ))
    table1Count.write("Gain/Loss\tFemale\tMale and Female\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)&(orthoExpressXA["xsome"]=="X")]),
        len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)&(orthoExpressXA["xsome"]=="A")]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_mel"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_mel"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_mel"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_mel"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_sim"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_sim"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_sim"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_sim"]==1)&(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_F_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
    ))
    table1Count.write("Gain/Loss\tMale and Female\tMale\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)&(orthoExpressXA["xsome"]=="X")]),
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)&(orthoExpressXA["xsome"]=="A")]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_M_sim"]==1)]),
    ))
    table1Count.write("Gain/Loss\tMale and Female\tFemale\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)&(orthoExpressXA["xsome"]=="X")]),
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)&(orthoExpressXA["xsome"]=="A")]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_F_sim"]==1)]),
    ))
    table1Count.write("Gain/Loss\tMale\tUnbiased\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_M_mel_U_sim"].sum() / orthoExpressXA["flag_M_mel_U_sim"].sum(),
    ))
    table1Count.write("Gain/Loss\tFemale\tUnbiased\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_F_mel_U_sim"].sum() / orthoExpressXA["flag_F_mel_U_sim"].sum(),
    ))
    table1Count.write("Gain/Loss\tMale and Female\tUnbiased\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)&(orthoExpressXA["xsome"]=="X")]),
        len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)&(orthoExpressXA["xsome"]=="A")]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_mel"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_sim"]==1)&(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_MF_mel"]==1)&(orthoExpressXA["flag_U_sim"]==1)]),
    ))
    table1Count.write("Gain/Loss\tUnbiased\tMale\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_M_sim_U_mel"].sum() / orthoExpressXA["flag_M_sim_U_mel"].sum(),
    ))
    table1Count.write("Gain/Loss\tUnbiased\tFemale\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="X"]["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["xsome"]=="A"]["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_mel"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_mel"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_mel"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_mel"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k4_sim"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k4_sim"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_male_k27_sim"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
        orthoExpressXA[orthoExpressXA["flag_has_female_k27_sim"]==1]["flag_F_sim_U_mel"].sum() / orthoExpressXA["flag_F_sim_U_mel"].sum(),
    ))
    table1Count.write("Gain/Loss\tUnbiased\tMale and Female\t{0} ({1}, {2})\t{3:.2%}\t{4:.2%}\t{5:.2%}\t{6:.2%}\t{7:.2%}\t{8:.2%}\t{9:.2%}\t{10:.2%}\n".format(
        len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)&(orthoExpressXA["xsome"]=="X")]),
        len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)&(orthoExpressXA["xsome"]=="A")]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_mel"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_mel"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_mel"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_mel"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k4_sim"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k4_sim"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_male_k27_sim"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
        len(orthoExpressXA[(orthoExpressXA["flag_has_female_k27_sim"]==1)&(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]) / len(orthoExpressXA[(orthoExpressXA["flag_U_mel"]==1)&(orthoExpressXA["flag_MF_sim"]==1)]),
    ))
    table1Count.write("Expressed\t{0} ({1}, {2})\t{3} ({4}, {5})\t{6} ({7}, {8})\t{9:.2%}\t{10:.2%}\t{11:.2%}\t{12:.2%}\t{13:.2%}\t{14:.2%}\t{15:.2%}\t{16:.2%}\n".format(
        len(melExpressXA),
        len(melExpressXA[melExpressXA["xsome"]=="X"]),
        len(melExpressXA[melExpressXA["xsome"]=="A"]),
        len(simExpressXA),
        len(simExpressXA[simExpressXA["xsome"]=="X"]),
        len(simExpressXA[simExpressXA["xsome"]=="A"]),
        len(orthoExpressXA),
        len(orthoExpressXA[orthoExpressXA["xsome"]=="X"]),
        len(orthoExpressXA[orthoExpressXA["xsome"]=="A"]),
        orthoExpressXA["flag_has_male_k4_mel"].sum() / len(orthoExpressXA),
        orthoExpressXA["flag_has_female_k4_mel"].sum() / len(orthoExpressXA),
        orthoExpressXA["flag_has_male_k27_mel"].sum() / len(orthoExpressXA),
        orthoExpressXA["flag_has_female_k27_mel"].sum() / len(orthoExpressXA),
        orthoExpressXA["flag_has_male_k4_sim"].sum() / len(orthoExpressXA),
        orthoExpressXA["flag_has_female_k4_sim"].sum() / len(orthoExpressXA),
        orthoExpressXA["flag_has_male_k27_sim"].sum() / len(orthoExpressXA),
        orthoExpressXA["flag_has_female_k27_sim"].sum() / len(orthoExpressXA),
    ))
    table1Count.close()

    # Get counts of chromatin presence for conserved/diverged genes
    expChromCount = open("{}/Table_1_counts.tsv".format(args.outDir),"w")
    expchrom = orthoExpressXA.copy()
    expchromConditions = [
        expchrom["flag_conserved_male_bias"] == 1,
        expchrom["flag_conserved_female_bias"] == 1,
        expchrom["flag_switch_M_F"] == 1,
        expchrom["flag_M_mel_U_sim"] == 1,
        expchrom["flag_F_mel_U_sim"] == 1,
        expchrom["flag_M_sim_U_mel"] == 1,
        expchrom["flag_F_sim_U_mel"] == 1,
    ]
    expchromChoices = [
        "Conserved Male",
        "Conserved Female",
        "Reversal",
        "Gain/Loss Mel Male",
        "Gain/Loss Mel Female",
        "Gain/Loss Sim Male",
        "Gain/Loss Sim Female",
    ]
    expchrom["exp"] = np.select(expchromConditions, expchromChoices, "Other")
    expchrom["chrom"] = np.where(
        expchrom["flag_has_male_k4_mel"] == 1,
        "X",
        "O"
    )
    expchrom["chrom"] = np.where(
        expchrom["flag_has_female_k4_mel"] == 1,
        expchrom["chrom"] + "X",
        expchrom["chrom"] + "O"
    )
    expchrom["chrom"] = np.where(
        expchrom["flag_has_male_k27_mel"] == 1,
        expchrom["chrom"] + "X",
        expchrom["chrom"] + "O"
    )
    expchrom["chrom"] = np.where(
        expchrom["flag_has_female_k27_mel"] == 1,
        expchrom["chrom"] + "X",
        expchrom["chrom"] + "O"
    )
    expchrom["chrom"] = np.where(
        expchrom["flag_has_male_k4_sim"] == 1,
        expchrom["chrom"] + "X",
        expchrom["chrom"] + "O"
    )
    expchrom["chrom"] = np.where(
        expchrom["flag_has_female_k4_sim"] == 1,
        expchrom["chrom"] + "X",
        expchrom["chrom"] + "O"
    )
    expchrom["chrom"] = np.where(
        expchrom["flag_has_male_k27_sim"] == 1,
        expchrom["chrom"] + "X",
        expchrom["chrom"] + "O"
    )
    expchrom["chrom"] = np.where(
        expchrom["flag_has_female_k27_sim"] == 1,
        expchrom["chrom"] + "X",
        expchrom["chrom"] + "O"
    )
    expchrom = expchrom.rename(columns={"chrom": "mMk4_mFk4_mMk27_mFk27_sMk4_sFk4_sMk27_sFk27"})
    pd.crosstab(expchrom["exp"], expchrom["mMk4_mFk4_mMk27_mFk27_sMk4_sFk4_sMk27_sFk27"])
    expChromCount.write("\t\n")
    expChromCount.write("X\t{},")
    expChromCount.close()

    if not args.inMapBias and not args.inConsMapBias:
        # Get counts of feature types detected by each chroamtin mark
        melFeatureCounts = pd.DataFrame(melFeature["featureType_x"].value_counts()).rename(columns={"featureType_x": "Total"})
        melFeatureCounts["D. melaongaster Male H3K4me3"] = melFeature.groupby("featureType_x")["flag_m_K4_on"].sum()
        melFeatureCounts["D. melaongaster Female H3K4me3"] = melFeature.groupby("featureType_x")["flag_f_K4_on"].sum()
        melFeatureCounts["D. melaongaster Male H3K27me2me3"] = melFeature.groupby("featureType_x")["flag_m_K27_on"].sum()
        melFeatureCounts["D. melaongaster Female H3K27me2me3"] = melFeature.groupby("featureType_x")["flag_f_K27_on"].sum()
        melFeaturePercent = melFeatureCounts[[c for c in melFeatureCounts.columns if c!="Total"]].div(melFeatureCounts["Total"], axis=0)
        simFeatureCounts = pd.DataFrame(simFeature["featureType_x"].value_counts()).rename(columns={"featureType_x": "Total"})
        simFeatureCounts["D. simulans Male H3K4me3"] = simFeature.groupby("featureType_x")["flag_m_K4_on"].sum()
        simFeatureCounts["D. simulans Female H3K4me3"] = simFeature.groupby("featureType_x")["flag_f_K4_on"].sum()
        simFeatureCounts["D. simulans Male H3K27me2me3"] = simFeature.groupby("featureType_x")["flag_m_K27_on"].sum()
        simFeatureCounts["D. simulans Female H3K27me2me3"] = simFeature.groupby("featureType_x")["flag_f_K27_on"].sum()
        simFeaturePercent = simFeatureCounts[[c for c in simFeatureCounts.columns if c!="Total"]].div(simFeatureCounts["Total"], axis=0)

        # Count sex-limited exonic features
        melFrag = melFeature[~(melFeature["featureID"].str.contains("3UTR")) & ~(melFeature["featureID"].str.contains("5UTR")) & ~(melFeature["featureID"].str.contains("intron")) & ~(melFeature["featureID"].str.contains("TSS")) & (~melFeature["featureID"].str.contains("intergenic"))].copy()
        simFrag = simFeature[~(simFeature["featureID"].str.contains("3UTR")) & ~(simFeature["featureID"].str.contains("5UTR")) & ~(simFeature["featureID"].str.contains("intron")) & ~(simFeature["featureID"].str.contains("TSS")) & (~simFeature["featureID"].str.contains("intergenic"))].copy()
        print("{0} out of {1} ({2:.2%}) D. melanogaster exonic features are sex-lmited ({3} in males, {4} in females)".format(
            len(melFrag[melFrag["flag_mel_m_on0_apn"]+melFrag["flag_mel_f_on0_apn"]==1]),
            len(melFrag),
            len(melFrag[melFrag["flag_mel_m_on0_apn"]+melFrag["flag_mel_f_on0_apn"]==1])/len(melFrag),
            len(melFrag[(melFrag["flag_mel_m_on0_apn"]==1)&(melFrag["flag_mel_f_on0_apn"]==0)]),
            len(melFrag[(melFrag["flag_mel_m_on0_apn"]==0)&(melFrag["flag_mel_f_on0_apn"]==1)])        
        ))
        print("{0} out of {1} ({2:.2%}) D. simulans exonic features are sex-lmited ({3} in males, {4} in females)".format(
            len(simFrag[simFrag["flag_sim_m_on0_apn"]+simFrag["flag_sim_f_on0_apn"]==1]),
            len(simFrag),
            len(simFrag[simFrag["flag_sim_m_on0_apn"]+simFrag["flag_sim_f_on0_apn"]==1])/len(simFrag),
            len(simFrag[(simFrag["flag_sim_m_on0_apn"]==1)&(simFrag["flag_sim_f_on0_apn"]==0)]),
            len(simFrag[(simFrag["flag_sim_m_on0_apn"]==0)&(simFrag["flag_sim_f_on0_apn"]==1)])        
        ))
    
        # Count sex-limited genes on all chromosomes and on just the X and autosomes
        print("{0} out of {1} ({2:.2%}) expressed genes sex-limited in D. melanogaster ({3} in males, {4} in females)".format(
            len(melExpressWsexLim[melExpressWsexLim["flag_sex_limited"]==1]),
            len(melExpressWsexLim),
            len(melExpressWsexLim[melExpressWsexLim["flag_sex_limited"]==1])/len(melExpressWsexLim),
            len(melExpressWsexLim[melExpressWsexLim["flag_male_limited"]==1]),
            len(melExpressWsexLim[melExpressWsexLim["flag_female_limited"]==1])
        ))
        melExpress2XA = melExpressWsexLim[melExpressWsexLim["xsome"].isin(["X", "A"])]
        print("Restricting to the X and autosomes only: {0} out of {1} ({2:.2%}) expressed genes sex-limited in D. melanogaster ({3} in males, {4} in females)".format(
            len(melExpress2XA[melExpress2XA["flag_sex_limited"]==1]),
            len(melExpress2XA),
            len(melExpress2XA[melExpress2XA["flag_sex_limited"]==1])/len(melExpress2XA),
            len(melExpress2XA[melExpress2XA["flag_male_limited"]==1]),
            len(melExpress2XA[melExpress2XA["flag_female_limited"]==1])
        ))
        print("{0} out of {1} ({2:.2%}) expressed genes sex-limited in D. simulans ({3} in males, {4} in females)".format(
            len(simExpressWsexLim[simExpressWsexLim["flag_sex_limited"]==1]),
            len(simExpressWsexLim),
            len(simExpressWsexLim[simExpressWsexLim["flag_sex_limited"]==1])/len(simExpressWsexLim),
            len(simExpressWsexLim[simExpressWsexLim["flag_male_limited"]==1]),
            len(simExpressWsexLim[simExpressWsexLim["flag_female_limited"]==1])
        ))
        simExpress2XA = simExpressWsexLim[simExpressWsexLim["xsome"].isin(["X", "A"])]
        print("Restricting to the X and autosomes only: {0} out of {1} ({2:.2%}) expressed genes sex-limited in D. simulans ({3} in males, {4} in females)".format(
            len(simExpress2XA[simExpress2XA["flag_sex_limited"]==1]),
            len(simExpress2XA),
            len(simExpress2XA[simExpress2XA["flag_sex_limited"]==1])/len(simExpress2XA),
            len(simExpress2XA[simExpress2XA["flag_male_limited"]==1]),
            len(simExpress2XA[simExpress2XA["flag_female_limited"]==1])
        ))
    #    del(simExpress2)
    #    del(simExpress2XA)    

###### ALL TESTS ARE BELOW

# Consistency of sex-biased expression with previous studies:
    # 1) One-to-one orthologous pairs with conserved male-bias and Newell 2016 male-biased input genes
    # 2) One-to-one orthologous pairs with conserved female-bias and Newell 2016 female-biased input genes
    # 3) One-to-one orthologous pairs with conserved male-bias and Newell 2016 male-biased TRAP genes
    # 4) One-to-one orthologous pairs with conserved female-bias and Newell 2016 female-biased TRAP genes
    # 5) One-to-one orthologous pairs with conserved male-bias and Garud 2015 genes in H12 peaks
    # 6) One-to-one orthologous pairs with conserved female-bias and Garud 2015 genes in H12 peaks
    # 7) One-to-one orthologous pairs with conserved male-biased and Innocenti and Morrow 2010 male fitness genes
    # 8) One-to-one orthologous pairs with conserved male-biased and Innocenti and Morrow 2010 antaonistic genes
    # 9) One-to-one orthologous pairs with conserved female-biased and Innocenti and Morrow 2010 female fitness genes
    # 10) One-to-one orthologous pairs with conserved female-biased and Innocenti and Morrow 2010 antaonistic genes

    # Output consistent with literature results
    consistantLitOut = open(
            "{}/expression_consistent_w_literature_tests.txt".format(args.outDir),
            "w"
    )

    consistantLitOut.write(
        "Consistency of sex-biased expression with previous studies:\n\n"
    )

    # Conserved sex-biased genes enriched for genes with alternative splicing (>1 transcript)
    ## !!! Add to supp file flag creation:
    orthoExpressXA["flag_multi_transcript_mel"] = np.where(orthoExpressXA["num_transcripts"]>1,1,0)
    orthoExpressXA["flag_multi_transcript_sim"] = np.where(orthoExpressXA["sim_num_transcripts"]>1,1,0)
    orthoExpressXA["flag_multi_transcript_both"] = np.where(orthoExpressXA["flag_multi_transcript_mel"]+orthoExpressXA["flag_multi_transcript_sim"]==2,1,0)
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_conserved_sex_bias",
        "flag_sex_det_pathway",
        consistantLitOut,
        "One-to-one orthologous pairs with conserved sex-bias and sex determination genes"
    )

    # Conserved sex-biased genes enriched for genes in the sex determination pathway
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_conserved_sex_bias",
        "flag_multi_transcript_both",
        consistantLitOut,
        "One-to-one orthologous pairs with conserved sex-bias and annotated >1 transcript in both species"
    )

# !!! Need to add to gene lists prior
#     # Get Mcintyre 2006 gene lists with dmel 617 ids
#     geneList = "~/mclab/SHARE/McIntyre_Lab/useful_dmel_data/gene_lists"
#     mcintyreS1 = pd.read_csv("{}/sex_bias_splice_McIntyre2006/mcintyre2006_Stb1_sig_sex_by_probe_617.csv".format(geneList))
#     mcintyreS5 = pd.read_csv("{}/sex_bias_splice_McIntyre2006/mcintyre2006_Stb5_sig_sex_by_probe_617.csv".format(geneList))
#     mcintyreDf = pd.merge(
#         mcintyreS1,
#         mcintyreS5,
#         how="outer",
#         on=["primary_FBgn", "symbol"],
#         validate="1:1"
#     )
#     # Merge in McIntyre 2006 list with mel gene-level sample flags
#     orthoExpressXA = pd.merge(
#         orthoExpressXA,
#         mcintyreDf,
#         how="outer",
#         left_on="mel_geneID",
#         right_on="primary_FBgn",
#         indicator="merge_check",
#         validate="m:1"
#     )
#     # Fix flag_mcintyre_sex_bias_2006
#     orthoExpressXA["flag_mcintyre2006_sex_bias"] = np.where(
#             (~orthoExpressXA["mcintyre2006_STb1_bias_towards"].isna())
#             | (~orthoExpressXA["mcintyre2006_STb5_bias_towards"].isna()),
#             1,
#             0
#     )
#     orthoExpressXA = orthoExpressXA[orthoExpressXA["merge_check"]!="right_only"].drop(columns=["merge_check"])

#     # Conserved sex-biased genes enriched for genes with sex-specific splicing in McIntyre 2006
#     test_2var(
#         orthoExpressXA.fillna(0),
#         "flag_conserved_sex_bias",
#         "flag_mcintyre2006_sex_bias",
#         consistantLitOut,
#         "One-to-one orthologous pairs with  with sex-specific splicing in McIntyre 2006"
#     )

    # Conserved male-biased genes enriched for Newell 2016 male-biased input genes
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_male_bias",
            "flag_Newell2016_m_bias_input",
            consistantLitOut,
            "One-to-one orthologous pairs with conserved male-bias and Newell 2016 male-biased input genes"
    )

    # Conserved female-biased genes enriched for Newell 2016 female-biased input genes
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_female_bias",
            "flag_Newell2016_f_bias_input",
            consistantLitOut,
            "One-to-one orthologous pairs with conserved female-bias and Newell 2016 female-biased input genes"
    )

    # Conserved male-biased genes enriched for Newell 2016 male-biased TRAP genes
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_male_bias",
            "flag_Newell2016_m_bias_TRAP",
            consistantLitOut,
            "One-to-one orthologous pairs with conserved male-bias and Newell 2016 male-biased TRAP genes"
    )

    # Conserved female-biased genes enriched for Newell 2016 female-biased trap genes
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_female_bias",
            "flag_Newell2016_f_bias_TRAP",
            consistantLitOut,
            "One-to-one orthologous pairs with conserved female-bias and Newell 2016 female-biased TRAP genes"
    )

    # Conserved male-biased genes enriched for Garud 2015 genes in H12 peaks
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_male_bias",
            "flag_garud_2015_h12_top50",
            consistantLitOut,
            "One-to-one orthologous pairs with conserved male-bias and Garud 2015 genes in H12 peaks"
    )

    # Conserved female-biased genes enriched for Garud 2015 genes in H12 peaks
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_female_bias",
            "flag_garud_2015_h12_top50",
            consistantLitOut,
            "One-to-one orthologous pairs with conserved female-bias and Garud 2015 genes in H12 peaks"
    )

    # Conserved unbiased genes enriched for Garud 2015 genes in H12 peaks
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_unbiased",
            "flag_garud_2015_h12_top50",
            consistantLitOut,
            "One-to-one orthologous pairs with conserved unbiased expression and Garud 2015 genes in H12 peaks"
    )

    # Divergent sex-biased genes enriched for Garud 2015 genes in H12 peaks
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_diverged_sex_bias",
            "flag_garud_2015_h12_top50",
            consistantLitOut,
            "One-to-one orthologous pairs with divergent sex-bias and Garud 2015 genes in H12 peaks"
    )

    # # Conserved male-biased genes enriched for Innocenti and Morrow 2010 male fitness genes
    # test_2var(
    #         orthoExpressXA.fillna(0),
    #         "flag_conserved_male_bias",
    #         "flag_innocenti_morrow_m_fitness",
    #         consistantLitOut,
    #         "One-to-one orthologous pairs with conserved male-bias and Innocenti and Morrow 2010 male fitness genes"
    # )

    # # Conserved male-biased genes enriched for Innocenti and Morrow 2010 antagonistic genes
    # test_2var(
    #         orthoExpressXA.fillna(0),
    #         "flag_conserved_male_bias",
    #         "flag_innocenti_morrow_antagonistic",
    #         consistantLitOut,
    #         "One-to-one orthologous pairs with conserved male-bias and Innocenti and Morrow 2010 antagonistic genes"
    # )

    # # Conserved female-biased genes enriched for Innocenti and Morrow 2010 female fitness genes
    # test_2var(
    #         orthoExpressXA.fillna(0),
    #         "flag_conserved_female_bias",
    #         "flag_innocenti_morrow_f_fitness",
    #         consistantLitOut,
    #         "One-to-one orthologous pairs with conserved female-bias and Innocenti and Morrow 2010 female fitness genes"
    # )

    # # Conserved female-biased genes enriched for Innocenti and Morrow 2010 antagonistic genes
    # test_2var(
    #         orthoExpressXA.fillna(0),
    #         "flag_conserved_female_bias",
    #         "flag_innocenti_morrow_antagonistic",
    #         consistantLitOut,
    #         "One-to-one orthologous pairs with conserved female-bias and Innocenti and Morrow 2010 antagonistic genes"
    # )

    # # Divergent sex-biased genes enriched for Innocenti and Morrow 2010 antagonistic genes
    # test_2var(
    #         orthoExpressXA.fillna(0),
    #         "flag_diverged_sex_bias",
    #         "flag_innocenti_morrow_antagonistic",
    #         consistantLitOut,
    #         "One-to-one orthologous pairs with divergent sex-bias and Innocenti and Morrow 2010 antagonistic genes"
    # )

    consistantLitOut.close()


# Consitency of conserved expression with fru regulation in previous studies:
    # 1) One-to-one orthologous pairs with conserved male-bias and Dalton 2013 genes regulated by fruM in males
    # 2) One-to-one orthologous pairs with conserved female-bias and Dalton 2013 genes regulated by fruM in males

    # Output consistent with literature results
    fruRegLitOut = open(
            "{}/fru_regulation_w_express_tests.txt".format(args.outDir),
            "w"
    )

    fruRegLitOut.write(
        "Consitency of conserved expression with previous studies:\n\n"
    )

    # 1-to-1 ortholog pair male-biased expression vs. Dalton 2013 fruM regulation in males
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_male_bias",
            "flag_dalton2013_fruM_male_regulation",
            fruRegLitOut,
            "One-to-one orthologous genes with male-biased expression vs. Dalton 2013 fruM regulation in males"
    )

    # 1-to-1 ortholog pair female-biased expression vs. Dalton 2013 fruM regulation in males
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_female_bias",
            "flag_dalton2013_fruM_male_regulation",
            fruRegLitOut,
            "One-to-one orthologous genes with female-biased expression vs. Dalton 2013 fruM regulation in males"
    )


    fruRegLitOut.close()

# Consitency of conserved expression with dsx regulation in previous studies:
    # 1) One-to-one orthologous pairs with conserved male-bias and Arbeitman 2016 genes regulated by dsx
    # 2) One-to-one orthologous pairs with conserved female-bias and Arbeitman 2016 genes regulated by dsx

    # Output consistent with literature results
    dsxRegLitOut = open(
            "{}/dsx_regulation_w_express_tests.txt".format(args.outDir),
            "w"
    )

    dsxRegLitOut.write(
        "Consitency of conserved expression with previous studies:\n\n"
    )

    # 1-to-1 ortholog pair male-biased expression vs. Arbeitman 2016 dsx regulation
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_male_bias",
            "flag_arbeitman2016_dsx_regulated",
            dsxRegLitOut,
            "One-to-one orthologous genes with male-biased expression vs. Arbeitman 2016 dsx regulation"
    )

    # 1-to-1 ortholog pair female-biased expression vs. Arbeitman 2016 dsx regulation
    test_2var(
            orthoExpressXA.fillna(0),
            "flag_conserved_female_bias",
            "flag_arbeitman2016_dsx_regulated",
            dsxRegLitOut,
            "One-to-one orthologous genes with female-biased expression vs. Arbeitman 2016 dsx regulation"
    )


    dsxRegLitOut.close()

# Enrichment of various results with positive selection in FlyDIVas

    # Output consistency with positive selection results
    posSelectOut = open(
        "{}/positive_selection_enrichment_tests.txt".format(args.outDir),
        'w'
    )

    posSelectOut.write(
        "Enrichment of various results with positive selection in FlyDIVas:\n\n"
    )

    posSelectDF = pd.DataFrame(columns=[flag for flag in orthoExpressXA.columns if "flag_flydivas" in  flag])
    posSelectDF.columns = pd.MultiIndex.from_product(
        [
            [flag.split("_")[2] for flag in orthoExpressXA.columns if "flag_flydivas" in  flag and "pos12" in flag],
            [flag.split("_")[3] for flag in orthoExpressXA.columns if "flag_flydivas" in  flag and "12spp" in flag]
        ], names=['Phylogeny', 'Model'])

    for flag in [c for c in orthoExpressXA.columns if "flag_flydivas" in c]:
        # 1-to-1 ortholog pair conserved male-biased expression vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_conserved_male_bias",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with conserved male-biased expression vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        if flag.split("_")[3] == "pos12":
            posSelectOut.write(
                "Of the {} conserved male-biased with significant {} pos12, "
                "{} are on the X and {} are on the autosomes:\n{}\n".format(
                    len(orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]==2]),
                    flag.split("_")[2],
                    len(orthoExpressXA[(orthoExpressXA["xsome"]=="X")&(orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]==2)]),
                    len(orthoExpressXA[(orthoExpressXA["xsome"]=="A")&(orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]==2)]),
                    orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]==2][["mel_geneSymbol", "xsome"]].to_string(index=False)
                )
            )
        posSelectDF.loc["Conserved Male-biased Expression",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]==2]
        )


        # 1-to-1 ortholog pair with conserved male-biased expression that have conserved presence of male open chromatin vs. each flyDivas flag
        test_2var(
                orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]==1].fillna(0),
                "flag_has_male_k4_both_species",
                flag,
                posSelectOut,
                "One-to-one orthologous conserved male-biased genes that have conserved presence of male open chromatin vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Male-biased Expression (Male H3K4me3 Both Species)",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]+orthoExpressXA["flag_has_male_k4_both_species"]==3]
        )

        # 1-to-1 ortholog pair with conserved male-biased expression that have presence of male open chromatin in mel only vs. each flyDivas flag
        test_2var(
                orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]==1].fillna(0),
                "flag_has_male_k4_mel_only",
                flag,
                posSelectOut,
                "One-to-one orthologous conserved male-biased genes that have presence of male open chromatin in D. melanogster only vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Male-biased Expression (Male H3K4me3 D. melanogaster only)",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]+orthoExpressXA["flag_has_male_k4_mel_only"]==3]
        )

                # 1-to-1 ortholog pair with conserved male-biased expression that have presence of male open chromatin in sim only vs. each flyDivas flag
        test_2var(
                orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]==1].fillna(0),
                "flag_has_male_k4_sim_only",
                flag,
                posSelectOut,
                "One-to-one orthologous conserved male-biased genes that have presence of male open chromatin in D. simulans only vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Male-biased Expression (Male H3K4me3 D. simulans only)",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA[flag]+orthoExpressXA["flag_has_male_k4_sim_only"]==3]
        )

        # 1-to-1 ortholog pair conserved female-biased expression vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_conserved_female_bias",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with conserved female-biased expression vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Female-biased Expression",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_conserved_female_bias"]+orthoExpressXA[flag]==2]
        )
        # 1-to-1 ortholog pair divergent sex-biased expression vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_diverged_sex_bias",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with divergent sex-biased expression vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Divergent Sex-biased Expression",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_diverged_sex_bias"]+orthoExpressXA[flag]==2]
        )
        # 1-to-1 ortholog pair mel-specific sex-biased expression vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_sex_bias_mel_only",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with D. melanogaster-specific sex-biased expression vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["D. melanogaster-specific Sex-biased Expression",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_sex_bias_mel_only"]+orthoExpressXA[flag]==2]
        )
        # 1-to-1 ortholog pair sim-specific sex-biased expression vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_sex_bias_sim_only",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with D. simulans-specific sex-biased expression vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["D. simulans-specific Sex-biased Expression",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_sex_bias_sim_only"]+orthoExpressXA[flag]==2]
        )
        # 1-to-1 ortholog pair female-biased in one species and male-biased in the other vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_switch_M_F",
                flag,
                posSelectOut,
                "One-to-one orthologous genes female-biased in one species and male-biased in the other vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Reversal of Sex-biased Expression",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_switch_M_F"]+orthoExpressXA[flag]==2]
        )
        if flag.split("_")[3] == "pos12":
            posSelectOut.write(
                "Of the {} divergence sex swithcing genes with significant {} pos12, "
                "{} are on the X and {} are on the autosomes:\n{}\n".format(
                    len(orthoExpressXA[orthoExpressXA["flag_switch_M_F"]+orthoExpressXA[flag]==2]),
                    flag.split("_")[2],
                    len(orthoExpressXA[(orthoExpressXA["xsome"]=="X")&(orthoExpressXA["flag_switch_M_F"]+orthoExpressXA[flag]==2)]),
                    len(orthoExpressXA[(orthoExpressXA["xsome"]=="A")&(orthoExpressXA["flag_switch_M_F"]+orthoExpressXA[flag]==2)]),
                    orthoExpressXA[orthoExpressXA["flag_switch_M_F"]+orthoExpressXA[flag]==2][["mel_geneSymbol", "xsome"]].to_string(index=False)
                )
            )
        # 1-to-1 ortholog pair female-biased in one species and unbiased in the other vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_switch_F_U",
                flag,
                posSelectOut,
                "One-to-one orthologous genes female-biased in one species and unbiased in the other vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Female-biased Expression in One Species",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_switch_F_U"]+orthoExpressXA[flag]==2]
        )
        # 1-to-1 ortholog pair male-biased in one species and unbiased in the other vs. each flyDivas flag
        test_2var(
                orthoExpressXA.fillna(0),
                "flag_switch_M_U",
                flag,
                posSelectOut,
                "One-to-one orthologous genes male-biased in one species and unbiased in the other vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Male-biased Expression in One Species",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoExpressXA[orthoExpressXA["flag_switch_M_U"]+orthoExpressXA[flag]==2]
        )

        # 1-to-1 ortholog pair conserved presence of male K4 vs. each flyDivas flag
        test_2var(
                orthoXA.fillna(0),
                "flag_has_male_k4_both_species",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with conserved presence of male H3K4me3 vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Presence of Male H3K4me3",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoXA[orthoXA["flag_has_male_k4_both_species"]+orthoXA[flag]==2]
        )
        # 1-to-1 ortholog pair conserved presence of female K4 vs. each flyDivas flag
        test_2var(
                orthoXA.fillna(0),
                "flag_has_female_k4_both_species",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with conserved presence of female H3K4me3 vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Presence of Female H3K4me3",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoXA[orthoXA["flag_has_female_k4_both_species"]+orthoXA[flag]==2]
        )
        # 1-to-1 ortholog pair conserved presence of male K27 vs. each flyDivas flag
        test_2var(
                orthoXA.fillna(0),
                "flag_has_male_k27_both_species",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with conserved presence of male H3K27me2me3 vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Presence of Male H3K27me2me3",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoXA[orthoXA["flag_has_male_k27_both_species"]+orthoXA[flag]==2]
        )
        # 1-to-1 ortholog pair conserved presence of female K27 vs. each flyDivas flag
        test_2var(
                orthoXA.fillna(0),
                "flag_has_female_k27_both_species",
                flag,
                posSelectOut,
                "One-to-one orthologous genes with conserved presence of female H3K27me2me3 vs. positive selection in {} {}".format(flag.split("_")[2], flag.split("_")[3])
        )
        posSelectDF.loc["Conserved Presence of Female H3K27me2me3",(
            flag.split("_")[2], flag.split("_")[3])] = len(
                orthoXA[orthoXA["flag_has_female_k27_both_species"]+orthoXA[flag]==2]
        )
        
    posSelectDF.to_csv("{}/positive_selection_counts.csv".format(args.outDir))
    posSelectOut.close()


# Chromatin and expression of sex-limited genes
    # Output chromatin and expression sex-limited comparison results
    sexLimOut = open(
        "{}/sex_limited_chromatin_expression_tests.txt".format(args.outDir),
        "w"
    )

    # D. melanogaster expressed genes with male-limited expression vs. female-limited expression
    test_2var(
            melExpressWsexLim[melExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited",
            "flag_female_limited",
            sexLimOut,
            "D. melanogaster genes with male-limited expression vs. female-limited expression on the X or autosomes."
    )

    # D. simulans expressed genes with male-limited expression vs. female-limited expression
    test_2var(
            simExpressWsexLim[simExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited",
            "flag_female_limited",
            sexLimOut,
            "D. simulans genes with male-limited expression vs. female-limited expression on the X or autosomes."
    )

    # In one-to-one orthologs with conserved male-limited expression vs.conserved female-limited expression
    test_2var(
            orthoExpressBothWsexLim[orthoExpressBothWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited_both_species",
            "flag_female_limited_both_species",
            sexLimOut,
            "One-to-one orthologs with conserved male-limited expression vs. conserved female-limited expression on the X or autosomes"
    )

    # D. melanogaster expressed genes with male-limited expression vs. male-limited H3K4me3
    test_2var(
            melExpressWsexLim[melExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited",
            "flag_male_limited_k4",
            sexLimOut,
            "D. melanogaster genes with male-limited expression vs. male-limited H3K4me3 on the X or autosomes"
    )
    # D. simulans expressed genes with male-limited expression vs. male-limited H3K4me3
    test_2var(
            simExpressWsexLim[simExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited",
            "flag_male_limited_k4",
            sexLimOut,
            "D. simulans genes with male-limited expression vs. male-limited H3K4me3 on the X or autosomes"
    )
    # In one-to-one orthologs with conserved male-limited expression vs.conserved male-limited H3K4me3
    test_2var(
            orthoExpressBothWsexLim[orthoExpressBothWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited_both_species",
            "flag_male_limited_k4_both_species",
            sexLimOut,
            "One-to-one orthologs with conserved male-limited expression vs. conserved male-limited H3K4me3 on the X or autosomes"
    )
    # D. melanogaster expressed genes with female-limited expression vs. female-limited H3K4me3
    test_2var(
            melExpressWsexLim[melExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_female_limited",
            "flag_female_limited_k4",
            sexLimOut,
            "D. melanogaster genes with female-limited expression vs. female-limited H3K4me3 on the X or autosomes"
    )
    # D. simulans expressed genes with female-limited expression vs. female-limited H3K4me3
    test_2var(
            simExpressWsexLim[simExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_female_limited",
            "flag_female_limited_k4",
            sexLimOut,
            "D. simulans genes with female-limited expression vs. female-limited H3K4me3 on the X or autosomes"
    )
    # In one-to-one orthologs with conserved female-limited expression vs.conserved female-limited H3K4me3
    test_2var(
            orthoExpressBothWsexLim[orthoExpressBothWsexLim["xsome"].isin(["X", "A"])],
            "flag_female_limited_both_species",
            "flag_female_limited_k4_both_species",
            sexLimOut,
            "One-to-one orthologs with conserved female-limited expression vs. conserved female-limited H3K4me3 on the X or autosomes"
    )

    # D. melanogaster expressed genes with male-limited expression vs. female-limited H3K27me2me3
    test_2var(
            melExpressWsexLim[melExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited",
            "flag_female_limited_k27",
            sexLimOut,
            "D. melanogaster genes with male-limited expression vs. female-limited H3K27me2me3 on the X or autosomes"
    )
    # D. simulans expressed genes with male-limited expression vs. female-limited H3K27me2me3
    test_2var(
            simExpressWsexLim[simExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited",
            "flag_female_limited_k27",
            sexLimOut,
            "D. simulans genes with male-limited expression vs. female-limited H3K27me2me3 on the X or autosomes"
    )
    # In one-to-one orthologs with conserved male-limited expression vs.conserved female-limited H3K27me2me3
    test_2var(
            orthoExpressBothWsexLim[orthoExpressBothWsexLim["xsome"].isin(["X", "A"])],
            "flag_male_limited_both_species",
            "flag_female_limited_k4_both_species",
            sexLimOut,
            "One-to-one orthologs with conserved male-limited expression vs. conserved female-limited H3K27me2me3 on the X or autosomes"
    )
    # D. melanogaster expressed genes with female-limited expression vs. male-limited H3K27me2me3
    test_2var(
            melExpressWsexLim[melExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_female_limited",
            "flag_male_limited_k27",
            sexLimOut,
            "D. melanogaster genes with female-limited expression vs. male-limited H3K27me2me3 on the X or autosomes"
    )
    # D. simulans expressed genes with female-limited expression vs. male-limited H3K27me2me3
    test_2var(
            simExpressWsexLim[simExpressWsexLim["xsome"].isin(["X", "A"])],
            "flag_female_limited",
            "flag_male_limited_k27",
            sexLimOut,
            "D. simulans genes with female-limited expression vs. male-limited H3K27me2me3 on the X or autosomes"
    )
    # In one-to-one orthologs with conserved female-limited expression vs.conserved male-limited H3K27me2me3
    test_2var(
            orthoExpressBothWsexLim[orthoExpressBothWsexLim["xsome"].isin(["X", "A"])],
            "flag_female_limited_both_species",
            "flag_male_limited_k27_both_species",
            sexLimOut,
            "One-to-one orthologs with conserved female-limited expression vs. conserved male-limited H3K27me2me3 on the X or autosomes"
    )

    sexLimOut.close()

# Chromatin and expression divergence tests
    # Output chromatin and expression sex bias results
    chromExpDivOut = open(
            "{}/chromatin_expression_divergence_tests.txt".format(args.outDir),
            "w"
    )

    chromExpDivOut.write(
        "Chromatin and expression divergence tests\n\n"
    )

    # 1-to-1 ortholog pairs with conserved male-biased expression and male H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]==1],
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with conserved male-biased expression and male H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with conserved male-biased expression and female H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]==1],
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with conserved male-biased expression and female H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with sim-specific male-biased expression and male H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_M_sim_U_mel"]==1],
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. simulans-specific male-biased expression and male H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with sim-specific male-biased expression and female H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_M_sim_U_mel"]==1],
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. simulans-specific male-biased expression and female H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with mel-specific male-biased expression and male H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_M_mel_U_sim"]==1],
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. melanogaster-specific male-biased expression and male H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with mel-specific male-biased expression and female H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_M_mel_U_sim"]==1],
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. melanogaster-specific male-biased expression and female H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with sim-specific female-biased expression and female H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_sim_U_mel"]==1],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. simulans-specific female-biased expression and female H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with sim-specific female-biased expression and male H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_sim_U_mel"]==1],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. simulans-specific female-biased expression and male H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with mel-specific female-biased expression and female H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_mel_U_sim"]==1],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. melanogaster-specific female-biased expression and female H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with mel-specific female-biased expression and male H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_mel_U_sim"]==1],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. melanogaster-specific female-biased expression and male H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with reversal of sex-biased expression and male H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_switch_M_F"]==1],
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with reversal of sex-biased expression and male H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with reversal of sex-biased expression and female H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_switch_M_F"]==1],
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with reversal of sex-biased expression and female H3K27me2me3 in D. melanogaster vs. D. simulans"
    )


    # 1-to-1 ortholog pairs with conserved female-biased expression and female H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_conserved_female_bias"]==1],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with conserved female-biased expression and female H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with conserved female-biased expression and male H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_conserved_female_bias"]==1],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with conserved female-biased expression and male H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with sim-specific female-biased expression and female H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_sim_U_mel"]==1],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. simulans-specific female-biased expression and female H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with sim-specific female-biased expression and male H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_sim_U_mel"]==1],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. simulans-specific female-biased expression and male H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with mel-specific female-biased expression and female H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_sim_U_mel"]==1],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. melanogaster-specific female-biased expression and female H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with mel-specific female-biased expression and male H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_F_mel_U_sim"]==1],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with D. melanogaster-specific female-biased expression and male H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with reversal of sex-biased expression and female H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_switch_M_F"]==1],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with reversal of sex-biased expression and female H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with reversal of sex-biased expression and male H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_switch_M_F"]==1],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with reversal of sex-biased expression and male H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog pairs with sex-biased expression in at least one species and male H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_sex_biased_mel"]+orthoExpressXA["flag_sex_biased_sim"]>0],
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with sex-biased expression in at least one species and male H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with sex-biased expression in at least one species and female H3K4me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_sex_biased_mel"]+orthoExpressXA["flag_sex_biased_sim"]>0],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with sex-biased expression in at least one species and female H3K4me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with sex-biased expression in at least one species and male H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_sex_biased_mel"]+orthoExpressXA["flag_sex_biased_sim"]>0],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with sex-biased expression in at least one species and male H3K27me2me3 in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog pairs with sex-biased expression in at least one species and female H3K27me2me3 in mel vs. sim
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_sex_biased_mel"]+orthoExpressXA["flag_sex_biased_sim"]>0],
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            chromExpDivOut,
            "One-to-one orthologous genes with sex-biased expression in at least one species and female H3K27me2me3 in D. melanogaster vs. D. simulans"
    )

    chromExpDivOut.close()

# Chromatin and expression association
    # 1) D. melanogaster gene expression vs. presence of H3K4me3
    # 2) D. simulans gene expression vs. presence of H3K4me3
    # 3) In one-to-one orthologs on the X, D. melanogaster female-biased expression vs. presence of H3K4me3 in females
    # 4) In one-to-one orthologs on the X, D. melanogaster female-biased expression vs. presence of H3K27me2me3 in males
    # 5) In one-to-one orthologs on the X, D. melanogaster male-biased expression vs. presence of H3K27me2me3 in females
    # 6) In one-to-one orthologs on the X, D. melanogaster male-biased expression vs. presence of H3K4me3 in males
    # 7) In one-to-one orthologs on the autosomes, D. melanogaster female-biased expression vs. presence of H3K4me3 in females
    # 8) In one-to-one orthologs on the autosomes, D. melanogaster female-biased expression vs. presence of H3K27me2me3 in males
    # 9) In one-to-one orthologs on the autosomes, D. melanogaster male-biased expression vs. presence of H3K27me2me3 in females
    # 10) In one-to-one orthologs on the autosomes, D. melanogaster male-biased expression vs. presence of H3K4me3 in males
    # 11) In one-to-one orthologs on the X, D. simulans female-biased expression vs. presence of H3K4me3 in females
    # 12) In one-to-one orthologs on the X, D. simulans female-biased expression vs. presence of H3K27me2me3 in males
    # 13) In one-to-one orthologs on the X, D. simulans male-biased expression vs. presence of H3K4me3 in females
    # 14) In one-to-one orthologs on the X, D. simulans male-biased expression vs. presence of H3K27me2me3 in females
    # 15) In one-to-one orthologs on the X, D. simulans male-biased expression vs. presence of H3K4me3 in males
    # 16) In one-to-one orthologs on the autosomes, D. simulans female-biased expression vs. presence of H3K4me3 in females
    # 17) In one-to-one orthologs on the autosomes, D. simulans female-biased expression vs. presence of H3K27me2me3 in males
    # 18) In one-to-one orthologs on the autosomes, D. simulans male-biased expression vs. presence of H3K27me2me3 in females
    # 19) In one-to-one orthologs on the autosomes, D. simulans male-biased expression vs. presence of H3K4me3 in males
    # 10) In one-to-one orthologs on the X or autosomes, D. melanogaster female-biased expression vs. presence of H3K4me3 in females
    # 21) In one-to-one orthologs on the X or autosomes, D. melanogaster female-biased expression vs. presence of H3K27me2me3 in males
    # 22) In one-to-one orthologs on the X or autosomes, D. melanogaster male-biased expression vs. presence of H3K27me2me3 in females
    # 23) In one-to-one orthologs on the X or autosomes, D. melanogaster male-biased expression vs. presence of H3K4me3 in males
    # 24) In one-to-one orthologs on the X or autosomes, D. simulans female-biased expression vs. presence of H3K4me3 in females
    # 25) In one-to-one orthologs on the X or autosomes, D. simulans female-biased expression vs. presence of H3K27me2me3 in males
    # 26) In one-to-one orthologs on the X or autosomes, D. simulans male-biased expression vs. presence of H3K27me2me3 in females
    # 27) In one-to-one orthologs on the X or autosomes, D. simulans male-biased expression vs. presence of H3K4me3 in males
    # 28) In one-to-one orthologs on the X or autosomes, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans
    # 29) In one-to-one orthologs on the X, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans
    # 30) In one-to-one orthologs on the autosomes, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans
    # 31) In one-to-one orthologs on the X or autosomes, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans
    # 32) In one-to-one orthologs on the X, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans
    # 33) In one-to-one orthologs on the autosomes, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans
    # 34) In one-to-one orthologs on the X or autosomes, male-biased expression with male K4 D. melanogaster vs. D. simulans
    # 35) In one-to-one orthologs on the X, male-biased expression with male K4 D. melanogaster vs. D. simulans
    # 36) In one-to-one orthologs on the autosomes, male-biased expression with male K4 D. melanogaster vs. D. simulans
    # 37) In one-to-one orthologs on the X or autosomes, male-biased expression with female K27 D. melanogaster vs. D. simulans
    # 38) In one-to-one orthologs on the X, male-biased expression with female K27 D. melanogaster vs. D. simulans
    # 39) In one-to-one orthologs on the autosomes, male-biased expression with female K27 D. melanogaster vs. D. simulans

    # Output chromatin sex bias results
    chromExpOut = open(
            "{}/chromatin_expression_tests.txt".format(args.outDir),
            "w"
    )

    chromExpOut.write(
        "Chromatin and expression association\n\n"
    )

    # Mel evidence of open chromatin in expressed genes
    test_2var(
            melXA,
            "flag_any_k4",
            "flag_expressed",
            chromExpOut,
            "D. melanogaster evidence of open chromatin in expressed genes"
    )

    # Mel evidence of open chromatin in expressed genes
    test_2var(
            simXA,
            "flag_any_k4",
            "flag_expressed",
            chromExpOut,
            "D. simulans evidence of open chromatin in expressed genes"
    )

    # Within all expressed genes or the one-to-one orthologs of each species, test expectation of chromatin in sex-biased expression
    testList = [
        "X female-biased female open chromatin",
        "A female-biased female open chromatin",
        "X female-biased male closed chromatin",
        "A female-biased male closed chromatin",
        "X male-biased male open chromatin",
        "A male-biased male open chromatin",
        "X male-biased female closed chromatin",
        "A male-biased female closed chromatin"        
    ]
    # Get genes that are one to one orthologs and conserved in sex bias or gain/loss sex bias
    melOrthoSBconservM = orthoDF[orthoDF["flag_conserved_male_bias"] + orthoDF["flag_conserved_unbiased"] > 0]["mel_geneID"].unique()
    simOrthoSBconservM = orthoDF[orthoDF["flag_conserved_male_bias"] + orthoDF["flag_conserved_unbiased"] > 0]["sim_geneID"].unique()
    melOrthoSBconservF = orthoDF[orthoDF["flag_conserved_female_bias"] + orthoDF["flag_conserved_unbiased"] > 0]["mel_geneID"].unique()
    simOrthoSBconservF = orthoDF[orthoDF["flag_conserved_female_bias"] + orthoDF["flag_conserved_unbiased"] > 0]["sim_geneID"].unique()
    melOrthoSBgainLossM = orthoDF[orthoDF["flag_M_mel_U_sim"] + orthoDF["flag_conserved_unbiased"] > 0]["mel_geneID"].unique()
    simOrthoSBgainLossM = orthoDF[orthoDF["flag_M_sim_U_mel"] + orthoDF["flag_conserved_unbiased"] > 0]["sim_geneID"].unique()
    melOrthoSBgainLossF = orthoDF[orthoDF["flag_F_mel_U_sim"] + orthoDF["flag_conserved_unbiased"] > 0]["mel_geneID"].unique()
    simOrthoSBgainLossF = orthoDF[orthoDF["flag_F_sim_U_mel"] + orthoDF["flag_conserved_unbiased"] > 0]["sim_geneID"].unique()

    for group in ["all", "orthologous", "conserved", "gainLoss"]:
        for species in ["melanogaster", "simulans"]:
            if species == "melanogaster":
                if group == "orthologous":
                    df = melExpress[melExpress["FBgn"].isin(melOrthoGenes)]
                else:
                    df = melExpress.copy()
            else:
                if group == "orthologous":
                    df = simExpress[simExpress["fbgn"].isin(simOrthoGenes)]
                else:
                    df = simExpress.copy()
            for chromosome in ["X", "A"]:
                if chromosome == "A":
                    name = "autosomes"
                else:
                    name = "X"
                for expression in ["female-biased", "male-biased"]:
                    if expression == "female-biased":
                        expflag = "flag_F"
                        if species == "melanogaster":
                            if group == "conserved":
                                chromExpress = df[df["FBgn"].isin(melOrthoSBconservF)]
                            elif group == "gainLoss":
                                chromExpress = df[df["FBgn"].isin(melOrthoSBgainLossF)]
                            else:
                                chromExpress = df.copy()
                        else:
                            if group == "conserved":
                                chromExpress = df[df["fbgn"].isin(simOrthoSBconservF)]
                            elif group == "gainLoss":
                                chromExpress = df[df["fbgn"].isin(simOrthoSBgainLossF)]  
                            else:
                                chromExpress = df.copy()
                    else:
                        expflag = "flag_M"
                        if species == "melanogaster":
                            if group == "conserved":
                                chromExpress = df[df["FBgn"].isin(melOrthoSBconservM)]
                            elif group == "gainLoss":
                                chromExpress = df[df["FBgn"].isin(melOrthoSBgainLossM)]
                            else:
                                chromExpress = df.copy()
                        else:
                            if group == "conserved":
                                chromExpress = df[df["fbgn"].isin(simOrthoSBconservM)]
                            elif group == "gainLoss":
                                chromExpress = df[df["fbgn"].isin(simOrthoSBgainLossM)]                            
                            else:
                                chromExpress = df.copy()
                    for chromatinSex in ["female", "male"]:
                        for chromatinType in ["open chromatin", "closed chromatin"]:
                            if chromatinType == "open chromatin":
                                chromflag = "flag_has_"+chromatinSex+"_k4"
                            else:
                                chromflag = "flag_has_"+chromatinSex+"_k27"
                            # Do tests in test list
                            if "{} {} {} {}".format(chromosome, expression, chromatinSex, chromatinType) in testList:
                                test_2var(
                                        chromExpress[chromExpress["xsome"]==chromosome],
                                        chromflag,
                                        expflag,
                                        chromExpOut,
                                        "In {} genes of D. {} presence of {} {} in {} genes on the {}".format(
                                            group,
                                            species,
                                            chromatinSex,
                                            chromatinType,
                                            expression,
                                            name
                                        )
                                )
                                crosstabTemp = pd.crosstab(chromExpress[chromExpress["xsome"]==chromosome][chromflag], chromExpress[chromExpress["xsome"]==chromosome][expflag])
                                chromExpOut.write("\t{0:.2%} genes with {1} expression have {2} {3}\n".format(
                                        crosstabTemp[1][1]/crosstabTemp[1].sum(),
                                        expression,
                                        chromatinSex,
                                        chromatinType
                                ))
                                chromExpOut.write("\t\t vs. {0:.2%} genes without {1} expression have {2} {3}\n\n".format(
                                        crosstabTemp[0][1]/crosstabTemp[0].sum(),
                                        expression,
                                        chromatinSex,
                                        chromatinType
                                ))

    #### For X and A combined
    for group in ["all", "orthologous", "conserved", "gainLoss"]:
        for species in ["melanogaster", "simulans"]:
            if species == "melanogaster":
                if group == "orthologous":
                    df = melExpressXA[melExpressXA["FBgn"].isin(melOrthoGenes)]
                else:
                    df = melExpressXA.copy()
            else:
                if group == "orthologous":
                    df = simExpressXA[simExpressXA["fbgn"].isin(simOrthoGenes)]
                else:
                    df = simExpressXA.copy()
            for expression in ["female-biased", "male-biased"]:
                if expression == "female-biased":
                    expflag = "flag_F"
                    if species == "melanogaster":
                        if group == "conserved":
                            chromExpress = df[df["FBgn"].isin(melOrthoSBconservF)]
                        elif group == "gainLoss":
                            chromExpress = df[df["FBgn"].isin(melOrthoSBgainLossF)]
                        else:
                            chromExpress = df.copy()
                    else:
                        if group == "conserved":
                            chromExpress = df[df["fbgn"].isin(simOrthoSBconservF)]
                        elif group == "gainLoss":
                            chromExpress = df[df["fbgn"].isin(simOrthoSBgainLossF)]  
                        else:
                            chromExpress = df.copy()
                else:
                    expflag = "flag_M"
                    if species == "melanogaster":
                        if group == "conserved":
                            chromExpress = df[df["FBgn"].isin(melOrthoSBconservM)]
                        elif group == "gainLoss":
                            chromExpress = df[df["FBgn"].isin(melOrthoSBgainLossM)]
                        else:
                            chromExpress = df.copy()
                    else:
                        if group == "conserved":
                            chromExpress = df[df["fbgn"].isin(simOrthoSBconservM)]
                        elif group == "gainLoss":
                            chromExpress = df[df["fbgn"].isin(simOrthoSBgainLossM)]                            
                        else:
                            chromExpress = df.copy()
                for chromatinSex in ["female", "male"]:
                    for chromatinType in ["open chromatin", "closed chromatin"]:
                        if chromatinType == "open chromatin":
                            chromflag = "flag_has_"+chromatinSex+"_k4"
                        else:
                            chromflag = "flag_has_"+chromatinSex+"_k27"
                        # Do tests in test list
                        if "{} {} {} {}".format(chromosome, expression, chromatinSex, chromatinType) in testList:
                            test_2var(
                                    chromExpress,
                                    chromflag,
                                    expflag,
                                    chromExpOut,
                                    "In {} genes of D. {} presence of {} {} in {} genes on the X or autosomes".format(
                                        group,
                                        species,
                                        chromatinSex,
                                        chromatinType,
                                        expression
                                    )
                            )
                            crosstabTemp = pd.crosstab(chromExpress[chromflag], chromExpress[expflag])
                            chromExpOut.write("\t{0:.2%} genes with {1} expression have {2} {3}\n".format(
                                    crosstabTemp[1][1]/crosstabTemp[1].sum(),
                                    expression,
                                    chromatinSex,
                                    chromatinType
                            ))
                            chromExpOut.write("\t\t vs. {0:.2%} genes without {1} expression have {2} {3}\n\n".format(
                                    crosstabTemp[0][1]/crosstabTemp[0].sum(),
                                    expression,
                                    chromatinSex,
                                    chromatinType
                            ))

    # In one-to-one orthologs on the X or autosomes, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA,
        "flag_M_OSCO_mel",
        "flag_M_OSCO_sim",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the X, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="X"],
        "flag_M_OSCO_mel",
        "flag_M_OSCO_sim",
        chromExpOut,
        "In one-to-one orthologs on the X, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the autosomes, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="A"],
        "flag_M_OSCO_mel",
        "flag_M_OSCO_sim",
        chromExpOut,
        "In one-to-one orthologs on the autosomes, male-biased expression OSCO (Mk4 or Fk27) D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the X or autosomes, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA,
        "flag_F_OSCO_mel",
        "flag_F_OSCO_sim",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the X, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="X"],
        "flag_F_OSCO_mel",
        "flag_F_OSCO_sim",
        chromExpOut,
        "In one-to-one orthologs on the X, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans"
    )
                
    # In one-to-one orthologs on the autosomes, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="A"],
        "flag_F_OSCO_mel",
        "flag_F_OSCO_sim",
        chromExpOut,
        "In one-to-one orthologs on the autosomes, female-biased expression OSCO (Fk4 or Mk27) D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the X or autosomes, male-biased expression with male K4 D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA,
        "flag_Mexp_Mk4_mel",
        "flag_Mexp_Mk4_sim",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, male-biased expression with male K4 D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the X, male-biased expression with male K4 D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="X"],
        "flag_Mexp_Mk4_mel",
        "flag_Mexp_Mk4_sim",
        chromExpOut,
        "In one-to-one orthologs on the X, male-biased expression with male K4 D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the autosomes, male-biased expression with male K4 D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="A"],
        "flag_Mexp_Mk4_mel",
        "flag_Mexp_Mk4_sim",
        chromExpOut,
        "In one-to-one orthologs on the autosomes, male-biased expression with male K4 D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the X or autosomes, male-biased expression with female K27 D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA,
        "flag_Mexp_Fk27_mel",
        "flag_Mexp_Fk27_sim",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, male-biased expression with female K27 D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the X, male-biased expression with female K27 D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="X"],
        "flag_Mexp_Fk27_mel",
        "flag_Mexp_Fk27_sim",
        chromExpOut,
        "In one-to-one orthologs on the X, male-biased expression with female K27 D. melanogaster vs. D. simulans"
    )

    # In one-to-one orthologs on the autosomes, male-biased expression with female K27 D. melanogaster vs. D. simulans
    test_2var(
        orthoExpressXA[orthoExpressXA["xsome"]=="X"],
        "flag_Mexp_Fk27_mel",
        "flag_Mexp_Fk27_sim",
        chromExpOut,
        "In one-to-one orthologs on the autosomes, male-biased expression with female K27 D. melanogaster vs. D. simulans"
    )

    # Are the 34 male-biased Fk27 genes enriched in the head of males or females?
    test_2var(
        orthoExpressXA[
            (orthoExpressXA["xsome"]=="X")
            & (orthoExpressXA["flag_flyAtlas2_Adult_Male_Head_expressed"]==1)].fillna(0),
        "flag_34_mbias_Fk27MelOnly_X",
        "flag_flyAtlas2_Adult_Male_Head_enrich_gt",
        chromExpOut,
        "In one-to-one orthologs on the X, conserved male-biased genes with female H3K27me2me3 only in D. melanogaster (34) vs. FlyAtlas2 Male Head enriched."
    )
    test_2var(
        orthoExpressXA[
            (orthoExpressXA["xsome"]=="X")
            & (orthoExpressXA["flag_flyAtlas2_Adult_Male_Head_expressed"]==1)].fillna(0),
        "flag_34_mbias_Fk27MelOnly_X",
        "flag_flyAtlas2_Adult_Female_Head_enrich_gt",
        chromExpOut,
        "In one-to-one orthologs on the X, conserved male-biased genes with female H3K27me2me3 only in D. melanogaster (34) vs. FlyAtlas2 Female Head enriched."
    )

    # In one-to-one orthologs on the X or autosomes with conserved unbiased expression vs. K4 in both sexes for of both species
    test_2var(
        orthoExpressXA,
        "flag_conserved_unbiased",
        "flag_both_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with conserved unbiased expression vs. K4 in both sexes for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with conserved unbiased expression vs. K4 in males for of both species
    test_2var(
        orthoExpressXA,
        "flag_conserved_unbiased",
        "flag_has_male_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with conserved unbiased expression vs. K4 in males for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with conserved unbiased expression vs. K4 in females for of both species
    test_2var(
        orthoExpressXA,
        "flag_conserved_unbiased",
        "flag_has_female_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with conserved unbiased expression vs. K4 in females for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with conserved male-biased expression vs. K4 in males for of both species
    test_2var(
        orthoExpressXA,
        "flag_conserved_male_bias",
        "flag_has_male_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with conserved male-biased expression vs. K4 in males for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with conserved female-biased expression vs. K4 in females for of both species
    test_2var(
        orthoExpressXA,
        "flag_conserved_female_bias",
        "flag_has_female_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with conserved female-biased expression vs. K4 in females for of both species"
    )

    # In one-to-one orthologs on the X or autosomes with conserved male-biased expression OR conserved unbiased vs. K4 in males for of both species
    test_2var(
        orthoExpressXA[orthoExpressXA["flag_conserved_male_bias"]+orthoExpressXA["flag_conserved_unbiased"]>0],
        "flag_conserved_male_bias",
        "flag_has_male_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with conserved male-biased expression OR conserved unbiased vs. K4 in males for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with conserved female-biased expression OR conserved unbiased vs. K4 in females for of both species
    test_2var(
        orthoExpressXA[orthoExpressXA["flag_conserved_female_bias"]+orthoExpressXA["flag_conserved_unbiased"]>0],
        "flag_conserved_female_bias",
        "flag_has_female_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with conserved female-biased expression OR conserved unbiased vs. K4 in females for of both species"
    )

    # In one-to-one orthologs on the X or autosomes with D. melanogaster-specific male-biased expression (sim unbiased) vs. K4 in males for of both species
    test_2var(
        orthoExpressXA,
        "flag_M_mel_U_sim",
        "flag_has_male_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with D. melanogaster-specific male-biased expression (sim unbiased) vs. K4 in males for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with D. simulans-specific male-biased expression (mel unbiased) vs. K4 in males for of both species
    test_2var(
        orthoExpressXA,
        "flag_M_sim_U_mel",
        "flag_has_male_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with D. simulans-specific male-biased expression (mel unbiased) vs. K4 in males for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with D. melanogaster-specific female-biased expression (sim unbiased) vs. K4 in fmales for of both species
    test_2var(
        orthoExpressXA,
        "flag_F_mel_U_sim",
        "flag_has_female_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with D. melanogaster-specific female-biased expression (sim unbiased) vs. K4 in females for of both species"
    )
    # In one-to-one orthologs on the X or autosomes with D. simulans-specific female-biased expression (mel unbiased) vs. K4 in females for of both species
    test_2var(
        orthoExpressXA,
        "flag_F_sim_U_mel",
        "flag_has_female_k4_both_species",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes with D. simulans-specific female-biased expression (mel unbiased) vs. K4 in females for of both species"
    )


    # In one-to-one orthologs on the X or autosomes, OSCO switch genes with Graze 2012 AI towards mel or sim
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_mel_MOSCO_sim_FOSCO",
        "flag_graze2012_AI_hybrid_any_mel_bias",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, D. melanogaster male OSCO and D. simuland female OSCO vs. Graze 2012 hybrid bias towards D. melanogaster alleles"
    )
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_mel_MOSCO_sim_FOSCO",
        "flag_graze2012_AI_hybrid_any_sim_bias",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, D. melanogaster male OSCO and D. simuland female OSCO vs. Graze 2012 hybrid bias towards D. simulans alleles"
    )
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_mel_FOSCO_sim_MOSCO",
        "flag_graze2012_AI_hybrid_any_mel_bias",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, D. melanogaster female OSCO and D. simuland male OSCO vs. Graze 2012 hybrid bias towards D. melanogaster alleles"
    )
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_mel_FOSCO_sim_MOSCO",
        "flag_graze2012_AI_hybrid_any_sim_bias",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, D. melanogaster female OSCO and D. simuland male OSCO vs. Graze 2012 hybrid bias towards D. simulans alleles"
    )
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_OSCO_switch",
        "flag_graze2012_AI_hybrid_any_mel_bias",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, OSCO switch (male in one species, female in the other) vs. Graze 2012 hybrid bias towards D. melanogaster alleles"
    )
    test_2var(
        orthoExpressXA.fillna(0),
        "flag_OSCO_switch",
        "flag_graze2012_AI_hybrid_any_sim_bias",
        chromExpOut,
        "In one-to-one orthologs on the X or autosomes, OSCO switch (male in one species, female in the other) vs. Graze 2012 hybrid bias towards D. simulans alleles"
    )




# Faster X-hypothesis: more sex-biased expression on X vs. A
#   1) One-to-one orthologous pairs with sex-bias in both species (regardless of direction) X vs. autosomes
#   2) One-to-one orthologous pairs with unbiased expression in both species X vs. autosomes
#   3) One-to-one orthologous pairs with conserved sex-bias (with conserved direction) X vs. autosomes
#   4) One-to-one orthologous pairs with conserved male-biased expression X vs. autosomes
#   5) One-to-one orthologous pairs with conserved female-biased expression X vs. autosomes

    # Output faster X results
    fasterXout = open(
            "{}/fasterX_X_A_tests.txt".format(args.outDir),
            "w"
    )

    fasterXout.write(
            "Faster X-hypothesis: more sex-biased expression on X vs. A\n\n"
    )

    # 1-to-1 ortholog sex-bias in both species X vs. A
    test_XA(
            orthoExpress,
            "flag_sex_bias_both_species",
            fasterXout,
            "One-to-one orthologous pairs with sex-bias in both species X vs. A"
    )

    # 1-to-1 ortholog unbiased in both species X vs. A
    test_XA(
            orthoExpress,
            "flag_conserved_unbiased",
            fasterXout,
            "One-to-one orthologous pairs with unbiased in both species X vs. A"
    )

    # 1-to-1 ortholog conserved sex-bias (conserved direction) X vs. A
    test_XA(
            orthoExpress,
            "flag_conserved_sex_bias",
            fasterXout,
            "One-to-one orthologous pairs with conserved sex-bias X vs. A"
    )

    # 1-to-1 ortholog conserved male-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_conserved_male_bias",
            fasterXout,
            "One-to-one orthologous pairs with conserved male-bias X vs. A"
    )

    # 1-to-1 ortholog conserved female-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_conserved_female_bias",
            fasterXout,
            "One-to-one orthologous pairs with conserved female-bias X vs. A"
    )

    # 1-to-1 ortholog mel-specific sex-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_sex_bias_mel_only",
            fasterXout,
            "One-to-one orthologous pairs with D. melanogaster-specific sex-bias X vs. A"
    )
    # 1-to-1 ortholog mel-specific male-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_M_mel_U_sim",
            fasterXout,
            "One-to-one orthologous pairs with D. melanogaster-specific male-bias X vs. A"
    )
    # 1-to-1 ortholog mel-specific female-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_F_mel_U_sim",
            fasterXout,
            "One-to-one orthologous pairs with D. melanogaster-specific female-bias X vs. A"
    )
    # 1-to-1 ortholog sim-specific sex-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_sex_bias_sim_only",
            fasterXout,
            "One-to-one orthologous pairs with D. simulans-specific sex-bias X vs. A"
    )
    # 1-to-1 ortholog sim-specific male-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_M_sim_U_mel",
            fasterXout,
            "One-to-one orthologous pairs with D. simulans-specific male-bias X vs. A"
    )
    # 1-to-1 ortholog sim-specific female-bias X vs. A
    test_XA(
            orthoExpress,
            "flag_F_sim_U_mel",
            fasterXout,
            "One-to-one orthologous pairs with D. simulans-specific female-bias X vs. A"
    )
    # 1-to-1 ortholog mel male-biased and sim female-biased X vs. A
    test_XA(
            orthoExpress,
            "flag_M_mel_F_sim",
            fasterXout,
            "One-to-one orthologous pairs with D. melaongaster male-bias and D. simulans female-bias X vs. A"
    )
    # 1-to-1 ortholog mel female-biased and sim male-biased X vs. A
    test_XA(
            orthoExpress,
            "flag_M_sim_F_mel",
            fasterXout,
            "One-to-one orthologous pairs with D. melaongaster female-bias and D. simulans male-bias X vs. A"
    )

    # In 1-to-1 orthologs Mel male-biased X vs. A
    test_XA(
            melExpress[melExpress["FBgn"].isin(melOrthoGenes)],
            "flag_M",
            fasterXout,
            "In orthologs of D. melanogaster male-biased X vs. autosomes"
    )
    
    # In 1-to-1 orthologs Mel female-biased X vs. A
    test_XA(
            melExpress[melExpress["FBgn"].isin(melOrthoGenes)],
            "flag_F",
            fasterXout,
            "In orthologs of D. melanogaster female-biased X vs. autosomes"
    )

    # In 1-to-1 orthologs Mel sex-biased X vs. A
    test_XA(
            melExpress[melExpress["FBgn"].isin(melOrthoGenes)],
            "flag_sex_biased",
            fasterXout,
            "In orthologs of D. melanogaster sex-biased X vs. autosomes"
    )

    # In 1-to-1 orthologs Sim male-biased X vs. A
    test_XA(
            simExpress[simExpress["fbgn"].isin(simOrthoGenes)],
            "flag_M",
            fasterXout,
            "In orthologs of D. simulans male-biased X vs. autosomes"
    )
    
    # In 1-to-1 orthologs Sim female-biased X vs. A
    test_XA(
            simExpress[simExpress["fbgn"].isin(simOrthoGenes)],
            "flag_F",
            fasterXout,
            "In orthologs of D. simulans female-biased X vs. autosomes"
    )

    # In 1-to-1 orthologs Sim sex-biased X vs. A
    test_XA(
            simExpress[simExpress["fbgn"].isin(simOrthoGenes)],
            "flag_sex_biased",
            fasterXout,
            "In orthologs of D. simulans sex-biased X vs. autosomes"
    )

    # Subset mel tests by flyAtlas2 head-enriched genes
    # In 1-to-1 orthologs Mel male-biased X vs. A
    test_XA(
            melExpress[
                (melExpress["FBgn"].isin(melOrthoGenes))
                & (melExpress["flag_flyAtlas2_Adult_Male_Head_enrich_gt"]==1)
            ],
            "flag_M",
            fasterXout,
            "In FlyAtlas2 male head-biased orthologs of D. melanogaster male-biased X vs. autosomes"
    )
    
    # In 1-to-1 orthologs Mel female-biased X vs. A
    test_XA(
            melExpress[
                (melExpress["FBgn"].isin(melOrthoGenes))
                & (melExpress["flag_flyAtlas2_Adult_Female_Head_enrich_gt"]==1)
            ],
            "flag_F",
            fasterXout,
            "In FlyAtlas2 female head-biased orthologs of D. melanogaster female-biased X vs. autosomes"
    )

    # In 1-to-1 orthologs Mel sex-biased X vs. A
    test_XA(
            melExpress[
                (melExpress["FBgn"].isin(melOrthoGenes))
                & (melExpress["flag_flyAtlas2_Adult_Female_Head_enrich_gt"] + melExpress["flag_flyAtlas2_Adult_Female_Head_enrich_gt"] > 0)
           ],
            "flag_sex_biased",
            fasterXout,
            "In FlyAtlas2 head-biased (either sex) orthologs of D. melanogaster sex-biased X vs. autosomes"
    )

    # Subset mel tests by flyAtlas2 non-head-enriched genes (head-enriched != 1 and head-expressed == 1)
    # In 1-to-1 orthologs Mel male-biased X vs. A
    test_XA(
            melExpress[
                (melExpress["FBgn"].isin(melOrthoGenes))
                & (melExpress["flag_flyAtlas2_Adult_Male_Head_enrich_gt"]!=1)
                & (melExpress["flag_flyAtlas2_Adult_Male_Head_expressed"]==1)
            ],
            "flag_M",
            fasterXout,
            "In FlyAtlas2 male nonhead-biased orthologs of D. melanogaster male-biased X vs. autosomes"
    )
    
    # In 1-to-1 orthologs Mel female-biased X vs. A
    test_XA(
            melExpress[
                (melExpress["FBgn"].isin(melOrthoGenes))
                & (melExpress["flag_flyAtlas2_Adult_Female_Head_enrich_gt"]!=1)
                & (melExpress["flag_flyAtlas2_Adult_Female_Head_expressed"]==1)
            ],
            "flag_F",
            fasterXout,
            "In FlyAtlas2 female non-head-biased orthologs of D. melanogaster female-biased X vs. autosomes"
    )

    # In 1-to-1 orthologs Mel sex-biased X vs. A
    test_XA(
            melExpress[
                (melExpress["FBgn"].isin(melOrthoGenes))
                & (melExpress["flag_flyAtlas2_Adult_Female_Head_enrich_gt"]!=1)
                & (melExpress["flag_flyAtlas2_Adult_Female_Head_enrich_gt"]!=1)
                & (
                    (melExpress["flag_flyAtlas2_Adult_Male_Head_expressed"]==1)
                    | (melExpress["flag_flyAtlas2_Adult_Female_Head_expressed"]==1)
                )
           ],
            "flag_sex_biased",
            fasterXout,
            "In FlyAtlas2 non-head-biased (either sex) orthologs of D. melanogaster sex-biased X vs. autosomes"
    )

    # In 1-to-1 orthologs conserved male-limited K4 X vs. A
    test_XA(
        orthoXA,
        "flag_male_limited_k4_both_species",
        fasterXout,
        "In one-to-one orthologs with conserved male-limited H3K4me3 X vs. autosomes"
    )
    # In 1-to-1 orthologs conserved female-limited K4 X vs. A
    test_XA(
        orthoXA,
        "flag_female_limited_k4_both_species",
        fasterXout,
        "In one-to-one orthologs with conserved female-limited H3K4me3 X vs. autosomes"
    )
    # In 1-to-1 orthologs conserved male-limited K27 X vs. A
    test_XA(
        orthoXA,
        "flag_male_limited_k27_both_species",
        fasterXout,
        "In one-to-one orthologs with conserved male-limited H3K27me2me3 X vs. autosomes"
    )
    # In 1-to-1 orthologs conserved female-limited K27 X vs. A
    test_XA(
        orthoXA,
        "flag_female_limited_k27_both_species",
        fasterXout,
        "In one-to-one orthologs with conserved female-limited H3K27me2me3 X vs. autosomes"
    )

    # In 1-to-1 orthologs of Mel any K4 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_any_k4",
        fasterXout,
        "In orthologs of D. melanogaster presence of any H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim any K4 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_any_k4",
        fasterXout,
        "In orthologs of D. simulans presence of any H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel any K27 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_any_k27",
        fasterXout,
        "In orthologs of D. melanogaster presence of any H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim any K27 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_any_k27",
        fasterXout,
        "In orthologs of D. simulans presence of any H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel male-limited K4 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_male_limited_k4",
        fasterXout,
        "In orthologs of D. melanogaster presence of male-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel female-limited K4 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_female_limited_k4",
        fasterXout,
        "In orthologs of D. melanogaster presence of female-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim male-limited K4 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_male_limited_k4",
        fasterXout,
        "In orthologs of D. simulans presence of male-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim female-limited K4 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_female_limited_k4",
        fasterXout,
        "In orthologs of D. simulans presence of female-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel male-limited K27 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_male_limited_k27",
        fasterXout,
        "In orthologs of D. melanogaster presence of male-limited H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel female-limited K27 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_female_limited_k27",
        fasterXout,
        "In orthologs of D. melanogaster presence of female-limited H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim male-limited K27 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_male_limited_k27",
        fasterXout,
        "In orthologs of D. simulans presence of male-limited H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim female-limited K27 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_female_limited_k27",
        fasterXout,
        "In orthologs of D. simulans presence of female-limited H3K27me2me3 X vs. A"
    )

    fasterXout.close()

## Expression sex bias of orthologs
#   1) One-to-one orthologous pairs with sex-bias on the X in D. melanogaster vs. D. simulans
#   2) One-to-one orthologous pairs with sex-bias on the autosomes in D. melanogaster vs. D. simulans
#   3) One-to-one orthologous pairs conserved male-biased vs. conserved female-biased on the X
#   4) One-to-one orthologous pairs conserved male-biased vs. conserved female-biased on the autosomes
#   5) One-to-one orthologous pairs conserved male-biased vs. conserved female-biased overall (X or A)
#   6) One-to-one orthologous pairs with unbiased D. melanogaster and sex-biased D. simulans, D. simulans-specific male-biased vs. female-biased overall (X or A)
#   7) One-to-one orthologous pairs with unbiased D. simulans and sex-biased D. melanogaster, D. melanogaster-specific male-biased vs. female-biased overall (X or A)
#   8) D. melanogaster male-biased vs. female-biased genes on the X
#   9) D. melanogaster male-biased vs. female-biased genes on the autosomes
#   10) D. simulans male-biased vs. female-biased genes on the X
#   11) D. simulans male-biased vs. female-biased genes on the autosomes
#   12) In orthologs of D. melanogaster male-biased vs. female-biased genes on the X
#   13) In orthologs of D. melanogaster male-biased vs. female-biased genes on the autosomes
#   14) In orthologs of D. simulans male-biased vs. female-biased genes on the X
#   15) In orthologs of D. simulans male-biased vs. female-biased genes on the autosomes

    # Output expression sex bias of orthologs results
    orthoExpOut = open(
            "{}/ortholog_expression_comparison_tests.txt".format(args.outDir),
            "w"
    )

    orthoExpOut.write(
            "Expression sex bias of orthologs\n\n"
    )

    # 1-to-1 ortholog detection of expression mel vs. sim
    test_2var(
        orthoExpressEitherXA,
        "flag_expressed_mel",
        "flag_expressed_sim",
        orthoExpOut,
        "One-to-one orthologous pairs with detection of expression on the X or autosomes in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog sex-bias on X mel vs. sim
    test_2var(
            orthoExpress[orthoExpress["xsome"]=="X"],
            "flag_sex_biased_mel",
            "flag_sex_biased_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with sex-bias on the X in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog sex-bias on autosomes mel vs. sim
    test_2var(
            orthoExpress[orthoExpress["xsome"]=="A"],
            "flag_sex_biased_mel",
            "flag_sex_biased_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with sex-bias on autosomes in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog male-bias on X or autosomes mel vs. sim
    test_2var(
            orthoExpressXA,
            "flag_M_mel",
            "flag_M_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with male-bias on X or autosomes in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog male-bias or unbiased on X or autosomes mel vs. sim
    test_2var(
            orthoExpressXA[
                (orthoExpressXA["flag_U_mel"]+orthoExpressXA["flag_M_mel"]==1)
                & (orthoExpressXA["flag_U_sim"]+orthoExpressXA["flag_M_sim"]==1)
            ],
            "flag_M_mel",
            "flag_M_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with male-bias or unbiased on X or autosomes in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog female-bias on X or autosomes mel vs. sim
    test_2var(
            orthoExpressXA,
            "flag_F_mel",
            "flag_F_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with female-bias on X or autosomes in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog female-bias or unbiased on X or autosomes mel vs. sim
    test_2var(
            orthoExpressXA[
                (orthoExpressXA["flag_U_mel"]+orthoExpressXA["flag_F_mel"]==1)
                & (orthoExpressXA["flag_U_sim"]+orthoExpressXA["flag_F_sim"]==1)
            ],
            "flag_F_mel",
            "flag_F_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with female-bias or unbiased on X or autosomes in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog mel-specific sex-biased on the X or autosomes male vs. female transitions
    test_2var(
            orthoExpressXA[
                (orthoExpressXA["flag_U_sim"]==1)
            ],
            "flag_M_mel",
            "flag_F_mel",
            orthoExpOut,
            "One-to-one orthologous pairs with D. melanogaster-specific sex-bias on X or autosomes in male vs. female transitions"
    )
    # 1-to-1 ortholog sim-specific sex-biased on the X or autosomes male vs. female transitions
    test_2var(
            orthoExpressXA[
                (orthoExpressXA["flag_U_mel"]==1)
            ],
            "flag_M_sim",
            "flag_F_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with D. simulans-specific sex-bias on X or autosomes in male vs. female transitions"
    )

    # 1-to-1 ortholog conserved male-biased vs. conserved female-biased on the X
    test_2var(
            orthoExpress[orthoExpress["xsome"]=="X"],
            "flag_conserved_male_bias",
            "flag_conserved_female_bias",
            orthoExpOut,
            "One-to-one orthologous pairs with conserved male-biased vs. conserved female-biased expression on the X"
    )

    # 1-to-1 ortholog conserved male-biased vs. conserved female-biased on the A
    test_2var(
            orthoExpress[orthoExpress["xsome"]=="A"],
            "flag_conserved_male_bias",
            "flag_conserved_female_bias",
            orthoExpOut,
            "One-to-one orthologous pairs with conserved male-biased vs. conserved female-biased expression on the autosomes"
    )

    # 1-to-1 ortholog conserved male-biased vs. conserved female-biased on X or A
    test_2var(
            orthoExpressXA,
            "flag_conserved_male_bias",
            "flag_conserved_female_bias",
            orthoExpOut,
            "One-to-one orthologous pairs with conserved male-biased vs. conserved female-biased expression on the X or autosomes"
    )

    # 1-to-1 ortholog with unbiased Mel and sex-biased Sim, Sim-specific male-biased vs. female-biased X or A
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_sex_bias_sim_only"]==1],
            "flag_M_sim",
            "flag_F_sim",
            orthoExpOut,
            "One-to-one orthologous pairs with unbiased D. melanogaster and sex-biased D. simulans, D. simulans-specific male-biased vs. female-biased expression on the X or autosomes"
    )

    # 1-to-1 ortholog with unbiased Sim and sex-biased Mel, Mel-specific male-biased vs. female-biased X or A
    test_2var(
            orthoExpressXA[orthoExpressXA["flag_sex_bias_mel_only"]==1],
            "flag_M_mel",
            "flag_F_mel",
            orthoExpOut,
            "One-to-one orthologous pairs with unbiased D. simulans and sex-biased D. melanogaster, D. melanogaster-specific male-biased vs. female-biased expression on the X or autosomes"
    )

    # Mel male-biased vs. female-biased on the X
    test_2var(
            melExpress[melExpress["xsome"]=="X"],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "D. melanogaster male-biased vs. female-biased expression on the X"
    )

    # Mel male-biased vs. female-biased on the A
    test_2var(
            melExpress[melExpress["xsome"]=="A"],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "D. melanogaster male-biased vs. female-biased expression on the autosomes"
    )

    # Sim male-biased vs. female-biased on the X
    test_2var(
            simExpress[simExpress["xsome"]=="X"],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "D. simulans male-biased vs. female-biased expression on the X"
    )

    # Sim male-biased vs. female-biased on the A
    test_2var(
            simExpress[simExpress["xsome"]=="A"],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "D. simulans male-biased vs. female-biased expression on the autosomes"
    )

    # In 1-to-1 orthologs Mel male-biased vs. female-biased on the X
    test_2var(
            melExpress[(melExpress["xsome"]=="X") & (melExpress["FBgn"].isin(melOrthoGenes))],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "In orthologs of D. melanogaster male-biased vs. female-biased expression on the X"
    )

    # In 1-to-1 orthologs Mel male-biased vs. female-biased on the A
    test_2var(
            melExpress[(melExpress["xsome"]=="A") & (melExpress["FBgn"].isin(melOrthoGenes))],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "In orthologs of D. melanogaster male-biased vs. female-biased expression on the autosomes"
    )

    # In 1-to-1 orthologs Sim male-biased vs. female-biased on the X
    test_2var(
            simExpress[(simExpress["xsome"]=="X") & (simExpress["fbgn"].isin(simOrthoGenes))],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "In orthologs of D. simulans male-biased vs. female-biased expression on the X"
    )

    # In 1-to-1 orthologs Sim male-biased vs. female-biased on the A
    test_2var(
            simExpress[(simExpress["xsome"]=="A") & (simExpress["fbgn"].isin(simOrthoGenes))],
            "flag_M",
            "flag_F",
            orthoExpOut,
            "In orthologs of D. simulans male-biased vs. female-biased expression on the autosomes"
    )


    orthoExpOut.close()


## Chromatin of orthologs
#   1) In one-to-one orthologs on the X, D. melanogaster presence of male H3K4me3 vs. presence of female H3K4me3
#   2) In one-to-one orthologs on the atuosomes, D. melanogaster presence of male H3K4me3 vs. presence of female H3K4me3
#   3) In one-to-one orthologs on the X, D. simulans presence of male H3K4me3 vs. presence of female H3K4me3
#   4) In one-to-one orthologs on the atuosomes, D. simulans presence of male H3K4me3 vs. presence of female H3K4me3
#   5) In one-to-one orthologs on the X, D. melanogaster presence of male H3K27me2me3 vs. presence of female H3K27me2me3
#   6) In one-to-one orthologs on the atuosomes, D. melanogaster presence of male H3K27me2me3 vs. presence of female H3K27me2me3
#   7) In one-to-one orthologs on the X, D. simulans presence of male H3K27me2me3 vs. presence of female H3K27me2me3
#   8) In one-to-one orthologs on the atuosomes, D. simulans presence of male H3K27me2me3 vs. presence of female H3K27me2me3
#   9) In one-to-one orthologs of D. melanogaster presence of H3K4me3 on X vs. autosomes
#   10) In one-to-one orthologs of D. simulans presence of H3K4me3 on X vs. autosomes
#   11) In one-to-one orthologs of D. melanogaster presence of H3K27me2me3 on X vs. autosomes
#   12) In one-to-one orthologs of D. simulans presence of H3K27me2me3 on X vs. autosomes
#   13) In one-to-one orthologs of D. melanogaster presence of male-limited H3K4me3 on X vs. autosomes
#   14) In one-to-one orthologs of D. melanogaster presence of female-limited H3K4me3 on X vs. autosomes
#   15) In one-to-one orthologs of D. simulans presence of male-limited H3K4me3 on X vs. autosomes
#   16) In one-to-one orthologs of D. simulans presence of female-limited H3K4me3 on X vs. autosomes
#   17) In one-to-one orthologs on the X, presence of male-limited H3K4me3 in D. melanogaster vs. D. simulans
#   18) In one-to-one orthologs on the autosomes, presence of male-limited H3K4me3 in D. melanogaster vs. D. simulans
#   19) In one-to-one orthologs on the X, presence of female-limited H3K4me3 in D. melanogaster vs. D. simulans
#   20) In one-to-one orthologs on the autosomes, presence of female-limited H3K4me3 in D. melanogaster vs. D. simulans
#   21) In one-to-one orthologs on the X, presence of male-limited H3K27me2me3 in D. melanogaster vs. D. simulans
#   22) In one-to-one orthologs on the autosomes, presence of male-limited H3K27me2me3 in D. melanogaster vs. D. simulans
#   23) In one-to-one orthologs on the X, presence of female-limited H3K27me2me3 in D. melanogaster vs. D. simulans
#   24) In one-to-one orthologs on the autosomes, presence of female-limited H3K27me2me3 in D. melanogaster vs. D. simulans
#   25) In one-to-one orthologs on X or autosomes with any open chromatin (either sex) in D. melanogaster vs D. simulans
#   26) In one-to-one orthologs on X or autosomes with any closed chromatin (either sex) in D. melanogaster vs D. simulans
#   27) In one-to-one orthologs on the, presence of male H3K4me3 in D. melanogaster vs. D. simulans

    # Output chromatin of orthologs results
    orthoChromOut = open(
            "{}/ortholog_chromatin_comparison_tests.txt".format(args.outDir),
            "w"
    )

    orthoChromOut.write(
        "Chromatin of orthologs\n\n"
    )

    # In 1-to-1 orthologs of Mel the presence of male K4 vs. female K4 on the X
    test_2var(
        melDF[(melDF["xsome"]=="X") & (melDF["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male H3K4me3 vs. presence of female H3K4me3 on the X"
    )
    # In 1-to-1 orthologs of Mel the presence of male K4 vs. female K4 on the autosomes
    test_2var(
        melDF[(melDF["xsome"]=="A") & (melDF["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male H3K4me3 vs. presence of female H3K4me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Mel the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        melXA[melXA["FBgn"].isin(melOrthoGenes)],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male H3K4me3 vs. presence of female H3K4me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Sim the presence of male K4 vs. female K4 on the X
    test_2var(
        simDF[(simDF["xsome"]=="X") & (simDF["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. simulans presence of male H3K4me3 vs. presence of female H3K4me3 on the X"
    )
    # In 1-to-1 orthologs of Sim the presence of male K4 vs. female K4 on the autosomes
    test_2var(
        simDF[(simDF["xsome"]=="A") & (simDF["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. simulans presence of male H3K4me3 vs. presence of female H3K4me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Sim the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        simXA[(simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. simulans presence of male H3K4me3 vs. presence of female H3K4me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Mel the presence of male K27 vs. female K27 on the X
    test_2var(
        melDF[(melDF["xsome"]=="X") & (melDF["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X"
    )
    # In 1-to-1 orthologs of Mel the presence of male K27 vs. female K27 on the autosomes
    test_2var(
        melDF[(melDF["xsome"]=="A") & (melDF["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Mel the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        melXA[(melXA["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Sim the presence of male K27 vs. female K27 on the X
    test_2var(
        simDF[(simDF["xsome"]=="X") & (simDF["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. simulans presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X"
    )
    # In 1-to-1 orthologs of Sim the presence of male K27 vs. female K27 on the autosomes
    test_2var(
        simDF[(simDF["xsome"]=="A") & (simDF["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. simulans presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Sim the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        simXA[(simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. simulans presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Mel with male-biased expression the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        melXA[(melXA["flag_M"]==1) & (melXA["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster with male-biased expression presence of male H3K4me3 vs. presence of female H3K4e3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Mel with male-biased expression the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        melXA[(melXA["flag_M"]==1) & (melXA["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster with male-biased expression presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Mel with female-biased expression the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        melXA[(melXA["flag_F"]==1) & (melXA["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster with female-biased expression presence of male H3K4me3 vs. presence of female H3K4e3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Mel with female-biased expression the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        melXA[(melXA["flag_F"]==1) & (melXA["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster with female-biased expression presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Mel with unbiased expression the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        melXA[(melXA["flag_U"]==1) & (melXA["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster with unbiased expression presence of male H3K4me3 vs. presence of female H3K4e3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Mel with male-biased expression the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        melXA[(melXA["flag_U"]==1) & (melXA["FBgn"].isin(melOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster with unbiased expression presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Sim with male-biased expression the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        simXA[(simXA["flag_M"]==1) & (simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. simulans with male-biased expression presence of male H3K4me3 vs. presence of female H3K4e3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Sim with male-biased expression the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        simXA[(simXA["flag_M"]==1) & (simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. simulans with male-biased expression presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Sim with female-biased expression the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        simXA[(simXA["flag_F"]==1) & (simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. simulans with female-biased expression presence of male H3K4me3 vs. presence of female H3K4e3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Sim with female-biased expression the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        simXA[(simXA["flag_F"]==1) & (simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. simulans with female-biased expression presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Sim with unbiased expression the presence of male K4 vs. female K4 on the X or autosomes
    test_2var(
        simXA[(simXA["flag_U"]==1) & (simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k4",
        "flag_has_female_k4",
        orthoChromOut,
        "In orthologs of D. simulans with unbiased expression presence of male H3K4me3 vs. presence of female H3K4e3 on the X or autosomes"
    )
    # In 1-to-1 orthologs of Sim with male-biased expression the presence of male K27 vs. female K27 on the X or autosomes
    test_2var(
        simXA[(simXA["flag_U"]==1) & (simXA["fbgn"].isin(simOrthoGenes))],
        "flag_has_male_k27",
        "flag_has_female_k27",
        orthoChromOut,
        "In orthologs of D. simulans with unbiased expression presence of male H3K27me2me3 vs. presence of female H3K27me2me3 on the X or autosomes"
    )


    # 1-to-1 ortholog on X or autosoomes with male-limited open chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_male_limited_k4_mel",
            "flag_male_limited_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes male-limited open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with male-limited open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_male_limited_k4_mel",
            "flag_male_limited_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X male-limited open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on A with male-limited open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_male_limited_k4_mel",
            "flag_male_limited_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on autosomes male-limited open chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or autosoomes with female-limited open chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_female_limited_k4_mel",
            "flag_female_limited_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes female-limited open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with female-limited open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_female_limited_k4_mel",
            "flag_female_limited_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X female-limited open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on A with female-limited open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_female_limited_k4_mel",
            "flag_female_limited_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on autosomes female-limited open chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or autosoomes with male-limited closed chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_male_limited_k27_mel",
            "flag_male_limited_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes male-limited closed chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with male-limited closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_male_limited_k27_mel",
            "flag_male_limited_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X male-limited closed chromatin in D. melanogaster vs. D. simulans"
    )
    
    # 1-to-1 ortholog on A with male-limited closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_male_limited_k27_mel",
            "flag_male_limited_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on autosomes male-limited closed chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or autosoomes with female-limited closed chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_female_limited_k27_mel",
            "flag_female_limited_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes female-limited closed chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with female-limited closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_female_limited_k27_mel",
            "flag_female_limited_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X female-limited closed chromatin in D. melanogaster vs. D. simulans"
    )
    
    # 1-to-1 ortholog on A with female-limited closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_female_limited_k27_mel",
            "flag_female_limited_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on autosomes female-limited closed chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or A with any open chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_any_k4_mel",
            "flag_any_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes any open chromatin in D. melanogaster vs. D. simulans"
    )
    
    # 1-to-1 ortholog on X or A with any closed chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_any_k27_mel",
            "flag_any_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes any closed chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or A with male open chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes male open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with male open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X male open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on A with male open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_has_male_k4_mel",
            "flag_has_male_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the autosomes male open chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or A with female open chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes female open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with female open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X female open chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on A with female open chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_has_female_k4_mel",
            "flag_has_female_k4_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the autosomes female open chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or A with male closed chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes male closed chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with male closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X male closed chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on A with male closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_has_male_k27_mel",
            "flag_has_male_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the autosomes male closed chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X or A with female closed chromatin in Mel vs Sim
    test_2var(
            orthoXA,
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X or autosomes female closed chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on X with female closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="X")],
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the X female closed chromatin in D. melanogaster vs. D. simulans"
    )
    # 1-to-1 ortholog on A with female closed chromatin in Mel vs Sim
    test_2var(
            orthoDF[(orthoDF["xsome"]=="A")],
            "flag_has_female_k27_mel",
            "flag_has_female_k27_sim",
            orthoChromOut,
            "One-to-one orthologous pairs on the autosomes female closed chromatin in D. melanogaster vs. D. simulans"
    )

    # 1-to-1 ortholog on X with conserved presence of sex-limited open chromatin in male vs female
    test_2var(
            orthoDF[(orthoDF["flag_sex_limited_k4_both_species"]==1)&(orthoDF["xsome"]=="X")],
            "flag_male_limited_k4_both_species",
            "flag_female_limited_k4_both_species",
            orthoChromOut,
            "One-to-one orthologous pairs on the X with conserved presence of sex-limited open chromatin in male vs female"
    )
    # 1-to-1 ortholog on A with conserved presence of sex-limited open chromatin in male vs female
    test_2var(
            orthoDF[(orthoDF["flag_sex_limited_k4_both_species"]==1)&(orthoDF["xsome"]=="A")],
            "flag_male_limited_k4_both_species",
            "flag_female_limited_k4_both_species",
            orthoChromOut,
            "One-to-one orthologous pairs on the autosomes with conserved presence of sex-limited open chromatin in male vs female"
    )

    # 1-to-1 ortholog on X with conserved presence of sex-limited closed chromatin in male vs female
    test_2var(
            orthoDF[(orthoDF["flag_sex_limited_k27_both_species"]==1)&(orthoDF["xsome"]=="X")],
            "flag_male_limited_k27_both_species",
            "flag_female_limited_k27_both_species",
            orthoChromOut,
            "One-to-one orthologous pairs on the X with conserved presence of sex-limited closed chromatin in male vs female"
    )
    # 1-to-1 ortholog on A with conserved presence of sex-limited closed chromatin in male vs female
    test_2var(
            orthoDF[(orthoDF["flag_sex_limited_k27_both_species"]==1)&(orthoDF["xsome"]=="A")],
            "flag_male_limited_k27_both_species",
            "flag_female_limited_k27_both_species",
            orthoChromOut,
            "One-to-one orthologous pairs on the autosomes with conserved presence of sex-limited closed chromatin in male vs female"
    )

    orthoChromOut.close()

## Chromatin K4 vs. K27 gene-level
#   1) In one-to-one orthologs on the X, D. melanogaster presence of any H3K4me3 vs. presence of any H3K27me2me3
#   2) In one-to-one orthologs on the autosomes, D. melanogaster presence of any H3K4me3 vs. presence of any H3K27me2me3
#   3) In one-to-one orthologs on the X, D. simulans presence of any H3K4me3 vs. presence of any H3K27me2me3
#   4) In one-to-one orthologs on the autosomes, D. simulans presence of any H3K4me3 vs. presence of any H3K27me2me3
#   5) In one-to-one orthologs on the X, D. melanogaster presence of sex-limited H3K4me3 vs. presence of sex-limited H3K27me2me3
#   6) In one-to-one orthologs on the autosomes, D. melanogaster presence of sex-limited H3K4me3 vs. presence of sex-limited H3K27me2me3
#   7) In one-to-one orthologs on the X, D. simulans presence of sex-limited H3K4me3 vs. presence of sex-limited H3K27me2me3
#   8) In one-to-one orthologs on the autosomes, D. simulans presence of sex-limited H3K4me3 vs. presence of sex-limited H3K27me2me3

    # Output chromatin K4 vs. K27 gene-level results
    chromK4vK27 = open(
            "{}/chromatin_K4_vs_K27_tests.txt".format(args.outDir),
            "w"
    )
    chromK4vK27.write(
        "Chromatin K4 vs. K27 gene-level\n\n"
    )

    # In 1-to-1 orthologs of Mel presence of male K4 vs. male K27 on the X
    test_2var(
            melDF[(melDF["xsome"]=="X") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_has_male_k4",
            "flag_has_male_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of male H3K4me3 vs. male H3K27me2me3 on the X"
    )
    # In 1-to-1 orthologs of Mel presence of male K4 vs. male K27 on the autosomes
    test_2var(
            melDF[(melDF["xsome"]=="A") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_has_male_k4",
            "flag_has_male_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of male H3K4me3 vs. male H3K27me2me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Mel presence of male K4 vs. male K27 on the autosomes
    test_2var(
            melXA[(melXA["FBgn"].isin(melOrthoGenes))],
            "flag_has_male_k4",
            "flag_has_male_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of male H3K4me3 vs. male H3K27me2me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Sim presence of male K4 vs. male K27 on the X
    test_2var(
            simDF[(simDF["xsome"]=="X") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_has_male_k4",
            "flag_has_male_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of male H3K4me3 vs. male H3K27me2me3 on the X"
    )
    # In 1-to-1 orthologs of Sim presence of male K4 vs. male K27 on the autosomes
    test_2var(
            simDF[(simDF["xsome"]=="A") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_has_male_k4",
            "flag_has_male_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of male H3K4me3 vs. male H3K27me2me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Sim presence of male K4 vs. male K27 on the autosomes
    test_2var(
            simXA[(simXA["fbgn"].isin(simOrthoGenes))],
            "flag_has_male_k4",
            "flag_has_male_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of male H3K4me3 vs. male H3K27me2me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Mel presence of female K4 vs. female K27 on the X
    test_2var(
            melDF[(melDF["xsome"]=="X") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_has_female_k4",
            "flag_has_female_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of female H3K4me3 vs. female H3K27me2me3 on the X"
    )
    # In 1-to-1 orthologs of Mel presence of female K4 vs. female K27 on the autosomes
    test_2var(
            melDF[(melDF["xsome"]=="A") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_has_female_k4",
            "flag_has_female_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of female H3K4me3 vs. female H3K27me2me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Mel presence of female K4 vs. female K27 on the autosomes
    test_2var(
            melXA[(melXA["FBgn"].isin(melOrthoGenes))],
            "flag_has_female_k4",
            "flag_has_female_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of female H3K4me3 vs. female H3K27me2me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Sim presence of female K4 vs. female K27 on the X
    test_2var(
            simDF[(simDF["xsome"]=="X") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_has_female_k4",
            "flag_has_female_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of female H3K4me3 vs. female H3K27me2me3 on the X"
    )
    # In 1-to-1 orthologs of Sim presence of female K4 vs. female K27 on the autosomes
    test_2var(
            simDF[(simDF["xsome"]=="A") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_has_female_k4",
            "flag_has_female_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of female H3K4me3 vs. female H3K27me2me3 on the autosomes"
    )
    # In 1-to-1 orthologs of Sim presence of female K4 vs. female K27 on the autosomes
    test_2var(
            simXA[(simXA["fbgn"].isin(simOrthoGenes))],
            "flag_has_female_k4",
            "flag_has_female_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of female H3K4me3 vs. female H3K27me2me3 on the X or autosomes"
    )

    # In 1-to-1 orthologs of Mel presence of any K4 vs. any K27 on the X
    test_2var(
            melDF[(melDF["xsome"]=="X") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_any_k4",
            "flag_any_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of any H3K4me3 vs. any H3K27me2me3 on the X"
    )
    # In 1-to-1 orthologs of Mel presence of any K4 vs. any K27 on the autosomes
    test_2var(
            melDF[(melDF["xsome"]=="A") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_any_k4",
            "flag_any_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of any H3K4me3 vs. any H3K27me2me3 on the autosomes"
    )

    # In 1-to-1 orthologs of Sim presence of any K4 vs. any K27 on the X
    test_2var(
            simDF[(simDF["xsome"]=="X") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_any_k4",
            "flag_any_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of any H3K4me3 vs. any H3K27me2me3 on the X"
    )

    # In 1-to-1 orthologs of Sim presence of any K4 vs. any K27 on the autosomes
    test_2var(
            simDF[(simDF["xsome"]=="X") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_any_k4",
            "flag_any_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of any H3K4me3 vs. any H3K27me2me3 on the autosomes"
    )

    # In 1-to-1 orthologs of Mel presence of sex-limited K4 vs. sex-limited K27 on the X
    test_2var(
            melDF[(melDF["xsome"]=="X") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_sex_limited_k4",
            "flag_sex_limited_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of sex-limited H3K4me3 vs. sex-limited H3K27me2me3 on the X"
    )

    # In 1-to-1 orthologs of Mel presence of sex-limited K4 vs. sex-limited K27 on the autosomes
    test_2var(
            melDF[(melDF["xsome"]=="A") & (melDF["FBgn"].isin(melOrthoGenes))],
            "flag_sex_limited_k4",
            "flag_sex_limited_k27",
            chromK4vK27,
            "In orthologs of D. melanogaster presence of sex-limited H3K4me3 vs. sex-limited H3K27me2me3 on the autosomes"
    )

    # In 1-to-1 orthologs of Sim presence of sex-limited K4 vs. sex-limited K27 on the X
    test_2var(
            simDF[(simDF["xsome"]=="X") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_sex_limited_k4",
            "flag_sex_limited_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of sex-limited H3K4me3 vs. sex-limited H3K27me2me3 on the X"
    )

    # In 1-to-1 orthologs of Sim presence of sex-limited K4 vs. sex-limited K27 on the autosomes
    test_2var(
            simDF[(simDF["xsome"]=="A") & (simDF["fbgn"].isin(simOrthoGenes))],
            "flag_sex_limited_k4",
            "flag_sex_limited_k27",
            chromK4vK27,
            "In orthologs of D. simulans presence of sex-limited H3K4me3 vs. sex-limited H3K27me2me3 on the autosomes"
    )

    chromK4vK27.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
