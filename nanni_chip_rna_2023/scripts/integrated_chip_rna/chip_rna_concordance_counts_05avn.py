#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get concordance of gene sex biased expression and chromatin from exonic regions")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="CSV file of gene-level annotation flags with combo flag gene_sex_bias_ttest_foldchange")
    parser.add_argument("-p", "--prefix", dest="prefix", required=True, help="Prefix for files with gene groups of interest")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output text file for counts")
    parser.add_argument("-d", "--output-directory", dest="outDir", required=True, help="Output directory for gene groups of interest")

    args = parser.parse_args()
    return args

def main():
    # Get input file
    annotDF = pd.read_csv(args.inFile, low_memory=False)

    # To ensure nan values are counted, convert to string
    annotDF['gene_ratio2_ttest'] = annotDF['gene_ratio2_ttest'].fillna("nan")
    annotDF['gene_trend_ttest'] = annotDF['gene_trend_ttest'].fillna("nan")
    
    # Open output file    
    of = open(args.outFile, 'w')
    
    # Values in gene_sex_bias_ttest_foldchange: male_and_female, male, female, unbiased, ttestLowFC, nan
    # Values in gene_k4_sex: male_and_female, male, fem, unb, none
    # Coordinates are [gene_sex_bias_ttest_foldchange, gene_k4_sex]
    
    # Get crosstab of expression to K4 sex bias
    trendRatioStack = pd.DataFrame(annotDF.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_k4_sex'].value_counts()).rename(columns={'gene_k4_sex':'count'})
    trendRatioWide = trendRatioStack.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    of.write("\nExpression vs. K4 sex bias (all genes, all chromosomes)\n{}\n".format(trendRatioWide))
    limitedDF = annotDF[(annotDF['gene_rna']=="fem")|(annotDF['gene_rna']=="male")]
    trendRatioStackLimited = pd.DataFrame(limitedDF.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_k4_sex'].value_counts()).rename(columns={'gene_k4_sex':'count'})
    trendRatioWideLimited = trendRatioStackLimited.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    of.write("\nExpression vs. K4 sex bias (sex-limited genes only, all chromosomes)\n{}\n".format(trendRatioWideLimited))
    limitedDFX = annotDF[((annotDF['gene_rna']=="fem")|(annotDF['gene_rna']=="male"))&(annotDF['xsome']=="X")]
    trendRatioStackLimitedX = pd.DataFrame(limitedDFX.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_k4_sex'].value_counts()).rename(columns={'gene_k4_sex':'count'})
    trendRatioWideLimitedX = trendRatioStackLimitedX.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    of.write("\nExpression vs. K4 sex bias (sex-limited genes only, X chromosome)\n{}\n".format(trendRatioWideLimitedX))
    limitedDFA = annotDF[((annotDF['gene_rna']=="fem")|(annotDF['gene_rna']=="male"))&(annotDF['xsome']=="A")]
    trendRatioStackLimitedA = pd.DataFrame(limitedDFA.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_k4_sex'].value_counts()).rename(columns={'gene_k4_sex':'count'})
    trendRatioWideLimitedA = trendRatioStackLimitedA.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    of.write("\nExpression vs. K4 sex bias (sex-limited genes only, Autosomes)\n{}\n".format(trendRatioWideLimitedA))
    biasDF = annotDF[(annotDF['gene_rna']!="fem")&(annotDF['gene_rna']!="male")&(annotDF['gene_rna']!="none")]
    trendRatioStackBiased = pd.DataFrame(biasDF.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_k4_sex'].value_counts()).rename(columns={'gene_k4_sex':'count'})
    try:
        trendRatioWideBiased = trendRatioStackBiased.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                                 columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    except:
        trendRatioWideBiased = trendRatioStackBiased.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                                 columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1)
    of.write("\nExpression vs. K4 sex bias (sex-biased genes only, sex-limited excluded, all chromosomes)\n{}\n".format(trendRatioWideBiased))
    biasDFX = annotDF[((annotDF['gene_rna']!="fem")&(annotDF['gene_rna']!="male")&(annotDF['gene_rna']!="none"))&(annotDF['xsome']=="X")]
    trendRatioStackBiasedX = pd.DataFrame(biasDFX.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_k4_sex'].value_counts()).rename(columns={'gene_k4_sex':'count'})
    try:
        trendRatioWideBiasedX = trendRatioStackBiasedX.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                                   columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    except:
        trendRatioWideBiasedX = trendRatioStackBiasedX.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                                   columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1)

    of.write("\nExpression vs. K4 sex bias (sex-biased genes only, sex-limited excluded, X chromosome)\n{}\n".format(trendRatioWideBiasedX))
    biasDFA = annotDF[((annotDF['gene_rna']!="fem")&(annotDF['gene_rna']!="male")&(annotDF['gene_rna']!="none"))&(annotDF['xsome']=="A")]
    trendRatioStackBiasedA = pd.DataFrame(biasDFA.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_k4_sex'].value_counts()).rename(columns={'gene_k4_sex':'count'})
    try:
        trendRatioWideBiasedA = trendRatioStackBiasedA.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                                   columns='gene_k4_sex', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    except:
        trendRatioWideBiasedA = trendRatioStackBiasedA.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                                   columns='gene_k4_sex', fill_value=0)
    of.write("\nExpression vs. K4 sex bias (sex-biased genes only, sex-limited excluded, Autosomes)\n{}\n".format(trendRatioWideBiasedA))
    
    
    # Get crosstab of expressed genes and xsome
    trendRatioXsomeStack = pd.DataFrame(annotDF.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['xsome'].value_counts()).rename(columns={'xsome':'count'})
    trendRatioXsomeWide = trendRatioXsomeStack.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                                                 columns='xsome', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
#    expressXsome = pd.crosstab(annotDF['gene_sex_bias_ttest_foldchange'],annotDF['xsome']).drop(['nan','ttestLowFC'])
 
    
    # Output the concordant male and female biased genes
    femConcord = annotDF[(annotDF['gene_trend_ttest']=='female')&(annotDF['gene_k4_sex']=='fem')]
    maleConcord = annotDF[(annotDF['gene_trend_ttest']=='male')&(annotDF['gene_k4_sex']=='male')]
    femX = femConcord['xsome'].value_counts().reindex(['X']).fillna(0).loc['X']
    femA = femConcord['xsome'].value_counts().reindex(['A']).fillna(0).loc['A']
    maleX = maleConcord['xsome'].value_counts().reindex(['X']).fillna(0).loc['X']
    maleA = maleConcord['xsome'].value_counts().reindex(['A']).fillna(0).loc['A']
    of.write("\n{} genes with evidence for female-based expression and concordant female K4: {}\n".format(
            len(femConcord),femConcord['symbol'].tolist()))
    of.write("\nChromosomal distribution of concordant female genes:\n\t{0} out of {1} expressed genes on X ({2:.5%})".format(
            femX,trendRatioXsomeWide['X'].sum(),femX/trendRatioXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on autosomes ({2:.5%})\n".format(
            femA,trendRatioXsomeWide['A'].sum(),femA/trendRatioXsomeWide['A'].sum()))
    of.write("\n{} genes with evidence for male-based expression and concordant male K4: {}\n".format(
            len(maleConcord),maleConcord['symbol'].tolist()))
    of.write("\nChromosomal distribution of concordant male genes:\n\t{0} out of {1} expressed genes on X ({2:.5%})".format(
            maleX,trendRatioXsomeWide['X'].sum(),maleX/trendRatioXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on autosomes ({2:.5%})\n".format(
            maleA,trendRatioXsomeWide['A'].sum(),maleA/trendRatioXsomeWide['A'].sum()))
    femConcord.to_csv(args.outDir+"/"+args.prefix+"_concordant_female.csv", index=False)
    maleConcord.to_csv(args.outDir+"/"+args.prefix+"_concordant_male.csv", index=False)

    # Output genes with male_and_female K4 and any sex-biased expression
    mfConcord = annotDF[(annotDF['gene_trend_ttest'].isin(["female","male","male_and_female"]))&
                        (annotDF['gene_k4_sex']=="male_and_female")]
    of.write("\n{} genes with male_and_female H3K4me3 and evidence for sex-biased expression: {}\n".format(len(mfConcord),mfConcord['symbol'].tolist()))
    mfConcord.to_csv(args.outDir+"/"+args.prefix+"_male_and_female_K4_sex_bias_express.csv", index=False)
    
    # 1) Concordant unbiased/sex-biased expression and K4 out of the total number of genes with detected K4
    totalK4 = len(annotDF[annotDF['gene_k4_sex']!="none"])
    of.write("\n1) Concordant unbiased/sex-biased expression and K4 out of the total number of genes with detected K4\n")
    of.write("{} genes with detected H3K4me3\n".format(totalK4))
    of.write("{0} ({1:.5%}) of genes with detected H3K4me3 have female-biased H3K4me3 and evidence for female-biased expression\n".format(
            trendRatioWide.loc[('female',slice(None)),'fem'].sum(),trendRatioWide.loc[('female',slice(None)),'fem'].sum()/totalK4))
    of.write("{0} ({1:.5%}) of genes with detected H3K4me3 have male-biased H3K4me3 and evidence for male-biased expression\n".format(
            trendRatioWide.loc[('male',slice(None)),'male'].sum(),trendRatioWide.loc[('male',slice(None)),'male'].sum()/totalK4))
    of.write("{0} ({1:.5%}) of genes with detected H3K4me3 have unbiased H3K4me3 and unbiased expression\n".format(
            trendRatioWide.loc[('unbiased',slice(None)),'unb'].sum(),trendRatioWide.loc[('unbiased',slice(None)),'unb'].sum()/totalK4))
    of.write("{0} ({1:.5%}) of genes with detected H3K4me3 have unbiased H3K4me3 and detected expression\n".format(
            trendRatioWide['unb'].sum()/totalK4,trendRatioWide['unb'].sum()))
    
    # 2) Concordant sex-biased expression and K4 out of the total number of sex-biased K4 genes
    #   (compare to Brown and Bachtrog 9% of sex-specific chromatin has sex-biased expression)
    totalK4sex = len(annotDF[~annotDF['gene_k4_sex'].isin(["none","unb"])])
    maleK4sex = len(annotDF[annotDF['gene_k4_sex']=="male"])
    femK4sex = len(annotDF[annotDF['gene_k4_sex']=="fem"])
    mfK4sex = len(annotDF[annotDF['gene_k4_sex']=="male_and_female"])
    of.write("\n2) Concordant evidence for sex-biased expression and K4 out of the total number of sex-biased K4 genes\n")
    of.write("\t(compare to Brown and Bachtrog 9% of sex-specific chromatin has sex-biased expression)\n")
    of.write("{} genes with evidence for sex-biased H3K4me3\n".format(totalK4sex))
    of.write("{0} ({1:.5%}) of genes with female-biased H3K4me3 have evidence for female-biased expression\n".format(
            trendRatioWide.loc[('female',slice(None)),'fem'].sum(),trendRatioWide.loc[('female',slice(None)),'fem'].sum()/femK4sex))
    of.write("{0} ({1:.5%}) of genes with sex-biased H3K4me3 have evidence for female-biased expression\n".format(
            trendRatioWide.loc[('female',slice(None)),'fem'].sum(),trendRatioWide.loc[('female',slice(None)),'fem'].sum()/totalK4sex))
    of.write("{0} ({1:.5%}) of genes with male-biased H3K4me3 have evidence for male-biased expression\n".format(
            trendRatioWide.loc[('male',slice(None)),'male'].sum(),trendRatioWide.loc[('male',slice(None)),'male'].sum()/maleK4sex))
    of.write("{0} ({1:.5%}) of genes with sex-biased H3K4me3 have evidence for male-biased expression\n".format(
            trendRatioWide.loc[('male',slice(None)),'male'].sum(),trendRatioWide.loc[('male',slice(None)),'male'].sum()/totalK4sex))
    tempValue = (trendRatioWide.loc[('female',slice(None)),'male_and_female'].sum()+trendRatioWide.loc[('male',slice(None)),
                                'male_and_female'].sum()+trendRatioWide.loc[('male_and_female',slice(None)),'male_and_female'].sum())
    of.write("{0} ({1:.5%}) of genes with male and female H3K4me3 within a gene have evidence for sex-biased expression\n".format(
            tempValue,tempValue/mfK4sex))
    
    # 3) Any expression (regardless of sex bias) and K4 out of the total number of expressed genes
    #   (compare to Negre 66% of expressed genes have K4 detected)
    totalRNA = len(annotDF[annotDF['gene_trend_ttest']!="nan"])
    of.write("\n3) Any expression (regardless of sex bias) and K4 out of the total number of expressed genes\n")
    of.write("\t(compare to Negre 66% of expressed genes have K4 detected)\n")
    of.write("{} expressed genes\n".format(totalRNA))
    of.write("{0} ({1:.5%}) of genes expressed genes have evidence for H3K4me3 (regardless of sex bias)\n".format(
            trendRatioWide.drop(columns=['none']).sum().sum(),
            trendRatioWide.drop(columns=['none']).sum().sum()/totalRNA))
    
    # 4) Concordant unbiased/sex-biased expression and K4 out of the total number of expressed genes
    of.write("\n4) Concordant unbiased/sex-biased expression and K4 out of the total number of expressed genes\n")
    of.write("{} expressed genes\n".format(totalRNA))
    of.write("{0} ({1:.5%}) of expressed genes have evidence for female-biased expression and female-biased H3K4me3\n".format(
            trendRatioWide.loc[('female',slice(None)),'fem'].sum(),trendRatioWide.loc[('female',slice(None)),'fem'].sum()/totalRNA))
    of.write("{0} ({1:.5%}) of expressed genes have evidence for male-biased expression and male-biased H3K4me3\n".format(
            trendRatioWide.loc[('male',slice(None)),'male'].sum(),trendRatioWide.loc[('male',slice(None)),'male'].sum()/totalRNA))
    of.write("{0} ({1:.5%}) of expressed genes have evidence for unbiased expression and unbiased H3K4me3\n".format(
            trendRatioWide.loc[('unbiased',slice(None)),'unb'].sum(),trendRatioWide.loc[('unbiased',slice(None)),'unb'].sum()/totalRNA))
    of.write("{0} ({1:.5%}) of expressed genes have evidence for unbiased expression and detected H3K4me3\n".format(
            trendRatioWide['unb'].sum(),trendRatioWide['unb'].sum()/totalRNA))
    
    # 5) Concordant unbiased/sex-biased expression and K4 out of total number of expressed genes with K4 detected
    totalK4rna = len(annotDF[(annotDF['gene_trend_ttest']!="nan")&(annotDF['gene_k4_sex']!="none")])
    of.write("\n5) Concordant unbiased/sex-biased expression and K4 out of total number of expressed genes with K4 detected\n")
    of.write("{} genes with detected H3K4me3 and expression\n".format(totalK4rna))
    of.write("{0} ({1:.5%}) of expressed genes with detected H3K4me3 have evidence for female-biased expression and female-biased H3K4me3\n".format(
            trendRatioWide.loc[('female',slice(None)),'fem'].sum(),trendRatioWide.loc[('female',slice(None)),'fem'].sum()/totalK4rna))
    of.write("{0} ({1:.5%}) of expressed genes with detected H3K4me3 have evidence for male-biased expression and male-biased H3K4me3\n".format(
            trendRatioWide.loc[('male',slice(None)),'male'].sum(),trendRatioWide.loc[('male',slice(None)),'male'].sum()/totalK4rna))
    of.write("{0} ({1:.5%}) of expressed genes with detected H3K4me3 have evidence for unbiased expression and unbiased H3K4me3\n".format(
            trendRatioWide.loc[('unbiased',slice(None)),'unb'].sum(),trendRatioWide.loc[('unbiased',slice(None)),'unb'].sum()/totalK4rna))
    of.write("{0} ({1:.5%}) of expressed genes with detected H3K4me3 have evidence for unbiased expression and detected H3K4me3\n".format(
            trendRatioWide['unb'].sum(),trendRatioWide['unb'].sum()/totalK4rna))
    
    of.close()
  
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()