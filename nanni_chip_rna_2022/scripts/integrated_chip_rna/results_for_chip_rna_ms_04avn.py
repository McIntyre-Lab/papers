#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from statsmodels.stats.contingency_tables import mcnemar
from statsmodels.stats.proportion import proportion_confint

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
    print_test_results(file, title, crosstabDF, chiSquare, cramersV, fisherExact2, fisherExactLess, fisherExactGreater, MN, kappa)

def print_test_results(file, title, crosstabDF, chiSquare, cramersV, fisherExact2, fisherExactLess, fisherExactGreater, MN, kappa = None):
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
        outText = outText + "\tKappa = {}\n\n\n".format(kappa)
    else:
        outText = outText + "\n\n"
    file.write(outText)


def main():
    # Get gene inputs
    melDF = pd.read_csv(args.inMG, low_memory=False)
    simDF = pd.read_csv(args.inSG, low_memory=False)
    orthoDF = pd.read_csv(args.inO, low_memory=False)

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
    orthoDF = orthoDF[orthoDF["flag_one2one_ortholog"]==1]
    melOrthoGenes = orthoDF["mel_geneID"].unique()
    simOrthoGenes = orthoDF["sim_geneID"].unique()
    
    # Drop ortholog pairs that the chromosome does not match (e.g. X in mel but A in sim)
    orthoDF = orthoDF[orthoDF["xsome"] == orthoDF["sim_xsome"]]

    # Drop genes without K4
    melK4 = melDF[melDF["flag_any_k4"]==1].copy()
    simK4 = simDF[simDF["flag_any_k4"]==1].copy()
    
    # Drop genes without K27
    melK27 = melDF[melDF["flag_any_k27"]==1].copy()
    simK27 = simDF[simDF["flag_any_k27"]==1].copy()

    # Drop genes not expressed
    melExpress2 = melDF[melDF["flag_expressed"] == 1].copy()
    simExpress2 = simDF[simDF["flag_expressed"] == 1].copy()
    orthoExpress2 = orthoDF[
            (orthoDF["flag_expressed_mel"]==1)
            & (orthoDF["flag_expressed_sim"]==1)
    ].copy()
    
    # Drop sex-limited genes
    melExpress = melExpress2[melExpress2["flag_sex_limited"] == 0]
    simExpress = simExpress2[simExpress2["flag_sex_limited"] == 0]
    orthoExpress = orthoExpress2[
            (orthoExpress2["flag_sex_limited_mel"]==0)
            & (orthoExpress2["flag_sex_limited_sim"]==0)
    ]

    # Get dataframe subsets of X and A (no chrom 4) only
    melXA = melDF[melDF["xsome"].isin(["X", "A"])]
    simXA = simDF[simDF["xsome"].isin(["X", "A"])]
    orthoXA = orthoDF[orthoDF["xsome"].isin(["X", "A"])]
    melExpressXA = melExpress[melExpress["xsome"].isin(["X", "A"])]
    simExpressXA = simExpress[simExpress["xsome"].isin(["X", "A"])]
    orthoExpressXA = orthoExpress[orthoExpress["xsome"].isin(["X", "A"])]

###### ALL TESTS ARE BELOW

# Consistency of sex-biased expression with previous studies:
    # 1) One-to-one orthologous pairs with conserved male-bias and Newell 2016 male-biased input genes
    # 2) One-to-one orthologous pairs with conserved female-bias and Newell 2016 female-biased input genes
    # 3) One-to-one orthologous pairs with conserved male-bias and Newell 2016 male-biased TRAP genes
    # 4) One-to-one orthologous pairs with conserved female-bias and Newell 2016 female-biased TRAP genes

    # Output consistent with literature results
    consistantLitOut = open(
            "{}/expression_consistent_w_literature_tests.txt".format(args.outDir),
            "w"
    )

    consistantLitOut.write(
        "Consistency of sex-biased expression with previous studies:\n\n"
    )

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

    # Within the one-to-one orthologs of each species, test expectation of chromatin in sex-biased expression
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
    for species in ["melanogaster", "simulans"]:
        if species == "melanogaster":
            df = melExpress[melExpress["FBgn"].isin(melOrthoGenes)]
        else:
            df = simExpress[simExpress["fbgn"].isin(simOrthoGenes)]
        for chromosome in ["X", "A"]:
            if chromosome == "A":
                name = "autosomes"
            else:
                name = "X"
            for expression in ["female-biased", "male-biased"]:
                if expression == "unbiased":
                    expflag = "flag_U"
                elif expression == "female-biased":
                    expflag = "flag_F"
                else:
                    expflag = "flag_M"
                for chromatinSex in ["female", "male"]:
                    for chromatinType in ["open chromatin", "closed chromatin"]:
                        if chromatinType == "open chromatin":
                            chromflag = "flag_has_"+chromatinSex+"_k4"
                        else:
                            chromflag = "flag_has_"+chromatinSex+"_k27"
                        # Do tests in test list
                        if "{} {} {} {}".format(chromosome, expression, chromatinSex, chromatinType) in testList:
                            test_2var(
                                    df[df["xsome"]==chromosome],
                                    chromflag,
                                    expflag,
                                    chromExpOut,
                                    "In orthologs of D. {} presence of {} {} in {} genes on the {}".format(
                                        species,
                                        chromatinSex,
                                        chromatinType,
                                        expression,
                                        name
                                    )
                            )
                            crosstabTemp = pd.crosstab(df[df["xsome"]==chromosome][chromflag], df[df["xsome"]==chromosome][expflag])
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
    for species in ["melanogaster", "simulans"]:
        if species == "melanogaster":
            df = melExpressXA[melExpressXA["FBgn"].isin(melOrthoGenes)]
        else:
            df = simExpressXA[simExpressXA["fbgn"].isin(simOrthoGenes)]
        for expression in ["female-biased", "male-biased"]:
            if expression == "unbiased":
                expflag = "flag_U"
            elif expression == "female-biased":
                expflag = "flag_F"
            else:
                expflag = "flag_M"
            for chromatinSex in ["female", "male"]:
                for chromatinType in ["open chromatin", "closed chromatin"]:
                    if chromatinType == "open chromatin":
                        chromflag = "flag_has_"+chromatinSex+"_k4"
                    else:
                        chromflag = "flag_has_"+chromatinSex+"_k27"
                    # Do tests in test list
                    if "{} {} {} {}".format(chromosome, expression, chromatinSex, chromatinType) in testList:
                        test_2var(
                                df,
                                chromflag,
                                expflag,
                                chromExpOut,
                                "In orthologs of D. {} presence of {} {} in {} genes on the X or autosomes".format(
                                    species,
                                    chromatinSex,
                                    chromatinType,
                                    expression
                                )
                        )
                        crosstabTemp = pd.crosstab(df[chromflag], df[expflag])
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

    chromExpOut.close()


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

    # In 1-to-1 orthologs of Mel any K4 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_any_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of any H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim any K4 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_any_k4",
        orthoChromOut,
        "In orthologs of D. simulans presence of any H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel any K27 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_any_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of any H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim any K27 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_any_k27",
        orthoChromOut,
        "In orthologs of D. simulans presence of any H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel male-limited K4 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_male_limited_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel female-limited K4 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_female_limited_k4",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of female-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim male-limited K4 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_male_limited_k4",
        orthoChromOut,
        "In orthologs of D. simulans presence of male-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim female-limited K4 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_female_limited_k4",
        orthoChromOut,
        "In orthologs of D. simulans presence of female-limited H3K4me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel male-limited K27 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_male_limited_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of male-limited H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Mel female-limited K27 X vs. A
    test_XA(
        melDF[melDF["FBgn"].isin(melOrthoGenes)],
        "flag_female_limited_k27",
        orthoChromOut,
        "In orthologs of D. melanogaster presence of female-limited H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim male-limited K27 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_male_limited_k27",
        orthoChromOut,
        "In orthologs of D. simulans presence of male-limited H3K27me2me3 X vs. A"
    )

    # In 1-to-1 orthologs of Sim female-limited K27 X vs. A
    test_XA(
        simDF[simDF["fbgn"].isin(simOrthoGenes)],
        "flag_female_limited_k27",
        orthoChromOut,
        "In orthologs of D. simulans presence of female-limited H3K27me2me3 X vs. A"
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
