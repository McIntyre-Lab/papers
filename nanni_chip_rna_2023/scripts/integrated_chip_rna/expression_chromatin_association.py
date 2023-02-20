#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get association of gene sex biased expression and the presence of chromatin marks")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="CSV file of gene-level annotation flags with combo flag gene_sex_bias_ttest_foldchange")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output text file for counts")

    args = parser.parse_args()
    return args

def main():
    # Get input file
    annotDF = pd.read_csv(args.inFile, low_memory=False)

    # Values in gene_sex_bias_ttest_foldchange: male_and_female, male, female, unbiased, ttestLowFC, nan
    # Values in gene_k4_sex: male_and_female, male, fem, unb, none
    # Coordinates are [gene_sex_bias_ttest_foldchange, gene_k4_sex]

    # To ensure nan values are counted, convert to string
    annotDF['gene_ratio2_ttest'] = annotDF['gene_ratio2_ttest'].fillna("nan")
    annotDF['gene_trend_ttest'] = annotDF['gene_trend_ttest'].fillna("nan")
    
    # Set flag for having a m/f k4 or k27 mark
    annotDF['flag_has_f_k4'] = np.where(annotDF['gene_k4_sex'].isin(['fem','male_and_female','unb']),1,0)
    annotDF['flag_has_m_k4'] = np.where(annotDF['gene_k4_sex'].isin(['male','male_and_female','unb']),1,0)
    annotDF['flag_has_f_k27'] = np.where(annotDF['gene_k27_sex'].isin(['fem','male_and_female','unb']),1,0)
    annotDF['flag_has_m_k27'] = np.where(annotDF['gene_k27_sex'].isin(['male','male_and_female','unb']),1,0)
    
    # Open output file    
    of = open(args.outFile, 'w')
    
    # Get sex-biased genes on the X
    biasDFX = annotDF[((annotDF['gene_rna']!="fem")&(annotDF['gene_rna']!="male")&(annotDF['gene_rna']!="none"))&(annotDF['xsome']=="X")].copy()

    # Flag overcompensated and undercompensated genes (using trend_ttest only)
    biasDFX['flag_overcompensated'] = np.where(biasDFX['gene_trend_ttest']=="male",1,0)
    biasDFX['flag_undercompensated'] = np.where(biasDFX['gene_trend_ttest']=="female",1,0)
    
    # Get crosstab of undercompensated and flag_has_f_k4
    undHasFk4 = pd.crosstab(biasDFX['flag_undercompensated'],biasDFX['flag_has_f_k4'])
    of.write("{}\n".format(undHasFk4))
    of.write("{0:.2%} of undercompensated genes have concordant female K4\n".format(undHasFk4.loc[1,1]/undHasFk4.loc[1,:].sum()))
    # Calculate Fisher exact test
    undHasFk4E = fisher_exact(undHasFk4)
    of.write("Odds Ratio = {}, Fisher exact pval = {}\n".format(undHasFk4E[0],undHasFk4E[1]))
    # Calculate chi2
#    undHasFk4Chi = chi2_contingency(undHasFk4)
    # Calculate Cramer's V
#    undHasFk4V = np.sqrt((undHasFk4Chi[0]/undHasFk4.sum().sum()))
#    of.write("Chi2 = {}, pval = {}, Cramer's V = {}\n".format(undHasFk4Chi[0],undHasFk4Chi[1],undHasFk4V))
    of.write("{0} of the {1} genes without concordant female k4 have concordant male k27\n\t({2:.2%} of the {3:.2%}, or {4:.2%}, of undercompensated genes without k4)\n\n".format(
        biasDFX.groupby(['flag_undercompensated','flag_has_f_k4'])['flag_has_m_k27'].sum().loc[1,0],
        undHasFk4.loc[1,0],biasDFX.groupby(['flag_undercompensated','flag_has_f_k4'])['flag_has_m_k27'].sum().loc[1,0]/undHasFk4.loc[1,0],
        undHasFk4.loc[1,0]/undHasFk4.loc[1,:].sum(),
        biasDFX.groupby(['flag_undercompensated','flag_has_f_k4'])['flag_has_m_k27'].sum().loc[1,0]/undHasFk4.loc[1,:].sum()))

    
    # Get crosstab of overcompensated and flag_has_m_k4
    overHasMk4 = pd.crosstab(biasDFX['flag_overcompensated'],biasDFX['flag_has_m_k4'])
    of.write("{}\n".format(overHasMk4))
    of.write("{0:.2%} of overcompensated genes have concordant male K4\n".format(overHasMk4.loc[1,1]/overHasMk4.loc[1,:].sum()))
    # Calculate Fisher exact test
    overHasMk4FE = fisher_exact(overHasMk4)
    of.write("Odds Ratio = {}, Fisher exact pval = {}\n".format(overHasMk4FE[0],overHasMk4FE[1]))
    # Calculate chi2
#    overHasMk4Chi = chi2_contingency(overHasMk4)
    # Calculate Cramer's V
#    overHasMk4V = np.sqrt((overHasMk4Chi[0]/overHasMk4.sum().sum()))
#    of.write("Chi2 = {}, pval = {}, Cramer's V = {}\n".format(overHasMk4Chi[0],overHasMk4Chi[1],overHasMk4V))
    of.write("{0} of the {1} genes without concordant male k4 have concordant female k27\n\t({2:.2%} of the {3:.2%}, or {4:.2%}, of overcompensated genes without k4)\n\n".format(
        biasDFX.groupby(['flag_overcompensated','flag_has_m_k4'])['flag_has_f_k27'].sum().loc[1,0],
        overHasMk4.loc[1,0],biasDFX.groupby(['flag_overcompensated','flag_has_m_k4'])['flag_has_f_k27'].sum().loc[1,0]/overHasMk4.loc[1,0],
        overHasMk4.loc[1,0]/overHasMk4.loc[1,:].sum(),
        biasDFX.groupby(['flag_overcompensated','flag_has_m_k4'])['flag_has_f_k27'].sum().loc[1,0]/overHasMk4.loc[1,:].sum()))



    # Get sex-biased genes on the Autosomes
    biasDFA = annotDF[((annotDF['gene_rna']!="fem")&(annotDF['gene_rna']!="male")&(annotDF['gene_rna']!="none"))&(annotDF['xsome']=="A")].copy()

    # Flag overcompensated and undercompensated genes (using trend_ttest only)
    biasDFA['flag_male_bias'] = np.where(biasDFA['gene_trend_ttest']=="male",1,0)
    biasDFA['flag_female_bias'] = np.where(biasDFA['gene_trend_ttest']=="female",1,0)
    
    # Get crosstab of female-biased and flag_has_f_k4
    femHasFk4 = pd.crosstab(biasDFA['flag_female_bias'],biasDFA['flag_has_f_k4'])
    of.write("{}\n".format(femHasFk4))
    of.write("{0:.2%} of female-biased genes have concordant female K4\n".format(femHasFk4.loc[1,1]/femHasFk4.loc[1,:].sum()))
    # Calculate Fisher exact test
    femHasFk4FE = fisher_exact(femHasFk4)
    of.write("Odds Ratio = {}, Fisher exact pval = {}\n".format(femHasFk4FE[0],femHasFk4FE[1]))
    # Calculate chi2
#    femHasFk4Chi = chi2_contingency(femHasFk4)
    # Calculate Cramer's V
#    femHasFk4V = np.sqrt((femHasFk4Chi[0]/femHasFk4.sum().sum()))
#    of.write("Chi2 = {}, pval = {}, Cramer's V = {}\n".format(femHasFk4Chi[0],femHasFk4Chi[1],femHasFk4V))
    of.write("{0} of the {1} genes without concordant female k4 have concordant male k27\n\t({2:.2%} of the {3:.2%}, or {4:.2%}, of female-biased genes without k4)\n\n".format(
        biasDFA.groupby(['flag_female_bias','flag_has_f_k4'])['flag_has_m_k27'].sum().loc[1,0],
        femHasFk4.loc[1,0],biasDFA.groupby(['flag_female_bias','flag_has_f_k4'])['flag_has_m_k27'].sum().loc[1,0]/femHasFk4.loc[1,0],
        femHasFk4.loc[1,0]/femHasFk4.loc[1,:].sum(),
        biasDFA.groupby(['flag_female_bias','flag_has_f_k4'])['flag_has_m_k27'].sum().loc[1,0]/femHasFk4.loc[1,:].sum()))

    
    # Get crosstab of overcompensated and flag_has_m_k4
    maleHasMk4 = pd.crosstab(biasDFA['flag_male_bias'],biasDFA['flag_has_m_k4'])
    of.write("{}\n".format(maleHasMk4))
    of.write("{0:.2%} of male-biased genes have concordant male K4\n".format(maleHasMk4.loc[1,1]/maleHasMk4.loc[1,:].sum()))
    # Calculate Fisher exact test
    maleHasMk4FE = fisher_exact(maleHasMk4)
    of.write("Odds Ratio = {}, Fisher exact pval = {}\n".format(maleHasMk4FE[0],maleHasMk4FE[1]))
    # Calculate chi2
#    maleHasMk4Chi = chi2_contingency(maleHasMk4)
    # Calculate Cramer's V
#    maleHasMk4V = np.sqrt((maleHasMk4Chi[0]/maleHasMk4.sum().sum()))
#    of.write("Chi2 = {}, pval = {}, Cramer's V = {}\n".format(maleHasMk4Chi[0],maleHasMk4Chi[1],maleHasMk4V))
    of.write("{0} of the {1} genes without concordant male k4 have concordant female k27\n\t({2:.2%} of the {3:.2%}, or {4:.2%}, of overcompensated genes without k4)\n\n".format(
        biasDFA.groupby(['flag_male_bias','flag_has_m_k4'])['flag_has_f_k27'].sum().loc[1,0],
        maleHasMk4.loc[1,0],biasDFA.groupby(['flag_male_bias','flag_has_m_k4'])['flag_has_f_k27'].sum().loc[1,0]/maleHasMk4.loc[1,0],
        maleHasMk4.loc[1,0]/maleHasMk4.loc[1,:].sum(),
        biasDFA.groupby(['flag_male_bias','flag_has_m_k4'])['flag_has_f_k27'].sum().loc[1,0]/maleHasMk4.loc[1,:].sum()))

    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()