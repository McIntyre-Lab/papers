#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Add single categorical flag for gene sex bias from ttest + foldchange flags and get counts of sex biased expression")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="CSV file of gene-level annotation flags")

    # Output data
    parser.add_argument("-c", "--counts", dest="outCounts", required=True, help="Output text file for counts of gene expression characterization")
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file with new categorical sex bias flag")

    args = parser.parse_args()
    return args

def main():
    # Get input file
    annotDF = pd.read_csv(args.inFile, low_memory=False)
    
    # Make flag for gene sex bias with ttest and >= 2 fold change
    #   based on gene_sex_bias assignment in gene_count_12amm.sas
    #   unbiased are required to have no significant ttest but fold change does not matter
    #   when ttest is significant but foldchange<1 then use trend variables
    annotDF['gene_ratio2_ttest'] = np.where(
            (annotDF['num_frags_ttest_1_ratio2_M']>0) & (annotDF['num_frags_ttest_1_ratio2_F']>0), "male_and_female",
            np.where((annotDF['num_frags_ttest_1_ratio2_M']>0) & (annotDF['num_frags_ttest_1_ratio2_F']==0), "male",
                     np.where((annotDF['num_frags_ttest_1_ratio2_M']==0) & (annotDF['num_frags_ttest_1_ratio2_F']>0), "female",
                              np.where((annotDF['num_frags_ttest_1_ratio2_M']==0) & (annotDF['num_frags_ttest_1_ratio2_F']==0)
                              & (annotDF['num_frags_ttest_1_ratio2_U']>0), "ttestLowFC",
                              np.where((annotDF['num_frags_ttest_1_ratio2_M']==0) & (annotDF['num_frags_ttest_1_ratio2_F']==0)
                              & (annotDF['num_frags_ttest_1_ratio2_U']==0) & (annotDF['num_frags_ttest_0_ratio2_M']+annotDF['num_frags_ttest_0_ratio2_F']+annotDF['num_frags_ttest_0_ratio2_U']>0),"unbiased",
                              np.where((annotDF['num_frags_ttest'].isna())|(annotDF['gene_sex_bias'].isna()),np.nan,"oops"))))))
    annotDF['gene_trend_ttest'] = np.where(
            (annotDF['num_frags_ttest_1_trend_M']>0) & (annotDF['num_frags_ttest_1_trend_F']>0), "male_and_female",
            np.where((annotDF['num_frags_ttest_1_trend_M']>0) & (annotDF['num_frags_ttest_1_trend_F']==0), "male",
                     np.where((annotDF['num_frags_ttest_1_trend_M']==0) & (annotDF['num_frags_ttest_1_trend_F']>0), "female",
                              np.where((annotDF['num_frags_ttest_1_trend_M']==0) & (annotDF['num_frags_ttest_1_trend_F']==0)
                              & (annotDF['num_frags_ttest_1_trend_U']>0), "ttestLowFC",
                              np.where((annotDF['num_frags_ttest_1_trend_M']==0) & (annotDF['num_frags_ttest_1_trend_F']==0)
                              & (annotDF['num_frags_ttest_1_trend_U']==0) & (annotDF['num_frags_ttest_0_trend_M']+annotDF['num_frags_ttest_0_trend_F']+annotDF['num_frags_ttest_0_trend_U']>0),"unbiased",
                              np.where((annotDF['num_frags_ttest'].isna())|(annotDF['gene_sex_bias'].isna()),np.nan,"oops"))))))
    annotDF.to_csv(args.outFile, index=False)

    # Gene sex bias counts
    of = open(args.outCounts, 'w')
    of.write("\ngene_ratio2_ttest\n{}\n".format(annotDF['gene_ratio2_ttest'].value_counts().to_string()))
    of.write("\ngene_trend_ttest\n{}\n".format(annotDF['gene_trend_ttest'].value_counts().to_string()))
    of.write("\n{}\n".format(pd.crosstab(annotDF['gene_ratio2_ttest'],annotDF['gene_trend_ttest'])))

    # Get crosstab of expressed genes (indexed by the trend and ratio2 pair) and xsome
    # Drop empty (nan) values from crosstab
    expressXsomeStack = pd.DataFrame(annotDF.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['xsome'].value_counts()).rename(columns={'xsome':'count'})
    expressXsomeWide = expressXsomeStack.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                     columns='xsome', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"))
    of.write("\n{}\n".format(expressXsomeWide.to_string()))
    
    # Get crosstab of gene sex bias (indexed by the trend and ratio2 pair) and rna detection
    # Drop empty (nan/none) from crosstab
    expressDetectStack = pd.DataFrame(annotDF.groupby(['gene_trend_ttest','gene_ratio2_ttest'])['gene_rna'].value_counts()).rename(columns={'gene_rna':'count'})
    expressDetectWide = expressDetectStack.pivot_table(index=['gene_trend_ttest','gene_ratio2_ttest'],
                     columns='gene_rna', fill_value=0).droplevel(0,axis=1).drop(index=("nan","nan"),columns=['none'])
    of.write("\n{}\n".format(expressDetectWide.to_string()))
    
    # Get crosstab of gene sex bias and detection (indexed by the trend, ratio2, detection tuple)
    #   and xsome - Drop empty (nan/none) from counts
    expressDetectXsomeStack = pd.DataFrame(annotDF.groupby(['gene_rna','gene_trend_ttest','gene_ratio2_ttest'])['xsome'].value_counts()).rename(columns={'xsome':'count'})
    expressDetectXsomeWide = expressDetectXsomeStack.pivot_table(index=['gene_rna','gene_trend_ttest','gene_ratio2_ttest'],
                     columns='xsome', fill_value=0).droplevel(0,axis=1).drop(index=('none','nan','nan'))
    of.write("\n{}\n".format(expressDetectXsomeWide.to_string()))

    # Output counts of the following by chromosomal location:
    # 1) evidence of sex-limited (detection and trend the same sex)
    # 2) pronounced sex-limited (detection, trend, and ratio2 all the same)
    # 3) evidence of sex-biased (detection is both, and trend is one sex)
    # 4) pronounced sex-biased (trend and ratio2 the same sex but detection is both)
    # 5) evidence of mixed sex-bias (detection is both and trend is male_and_female)
    # 6) pronounced mixed sex-bias (detection is both, trend and ratio2 are male_and_female)
    fLimitedXevidence = expressDetectXsomeWide.loc[('fem','female'),'X'].sum()
    fLimitedAevidence = expressDetectXsomeWide.loc[('fem','female'),'A'].sum()
    fLimitedXpronounced = expressDetectXsomeWide.loc[('fem','female','female'),'X'].sum()
    fLimitedApronounced = expressDetectXsomeWide.loc[('fem','female','female'),'A'].sum()
    fBiasedXevidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'female'),'X'].sum()
    fBiasedAevidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'female'),'A'].sum()
    fBiasedXpronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'female','female'),'X'].sum()
    fBiasedApronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'female','female'),'A'].sum()
    mLimitedXevidence = expressDetectXsomeWide.loc[('male','male'),'X'].sum()
    mLimitedAevidence = expressDetectXsomeWide.loc[('male','male'),'A'].sum()
    mLimitedXpronounced = expressDetectXsomeWide.loc[('male','male','male'),'X'].sum()
    mLimitedApronounced = expressDetectXsomeWide.loc[('male','male','male'),'A'].sum()
    mBiasedXevidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male'),'X'].sum()
    mBiasedAevidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male'),'A'].sum()
    mBiasedXpronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male','male'),'X'].sum()
    mBiasedApronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male','male'),'A'].sum()
    mfXevidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male_and_female'),'X'].sum()
    mfXpronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male_and_female','male_and_female'),'X'].sum()
    mfAevidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male_and_female'),'A'].sum()
    mfApronounced =expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male_and_female','male_and_female'),'A'].sum()
    totalFlimitedEvidence = expressDetectXsomeWide.loc[('fem','female'),:].sum().sum()
    totalFlimitedPronounced = expressDetectXsomeWide.loc[('fem','female','female'),:].sum().sum()
    totalMlimitedEvidence = expressDetectXsomeWide.loc[('male','male'),:].sum().sum()
    totalMlimitedPronounced = expressDetectXsomeWide.loc[('male','male','male'),:].sum().sum()
    totalFbiasedEvidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'female'),:].sum().sum()
    totalFbiasedPronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'female','female'),:].sum().sum()
    totalMbiasedEvidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male'),:].sum().sum()
    totalMbiasedPronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male','male'),:].sum().sum()
    totalMFevidence = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male_and_female'),:].sum().sum()
    totalMFpronounced = expressDetectXsomeWide[~expressDetectXsomeWide.index.get_level_values(
            'gene_rna').isin(['fem','male'])].loc[(slice(None),'male_and_female','male_and_female'),:].sum().sum()
    # Output female counts
    # 1) evidence of sex-limited (detection and trend the same sex)
    of.write("\n{0} genes with evidence of female-limited:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalFlimitedEvidence,fLimitedXevidence,expressXsomeWide['X'].sum(),
            fLimitedXevidence/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            fLimitedAevidence,expressXsomeWide['A'].sum(),fLimitedAevidence/expressXsomeWide['A'].sum()))
    # 2) pronounced sex-limited (detection, trend, and ratio2 all the same)
    of.write("\n{0} genes with pronounced female-limited:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalFlimitedPronounced,fLimitedXpronounced,expressXsomeWide['X'].sum(),
            fLimitedXpronounced/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            fLimitedApronounced,expressXsomeWide['A'].sum(),fLimitedApronounced/expressXsomeWide['A'].sum()))
    # 3) evidence of sex-biased (detection is both, and trend is one sex)
    of.write("\n{0} genes with evidence of female-bias:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalFbiasedEvidence,fBiasedXevidence,expressXsomeWide['X'].sum(),
            fBiasedXevidence/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            fBiasedAevidence,expressXsomeWide['A'].sum(),fBiasedAevidence/expressXsomeWide['A'].sum()))
    # 4) pronounced sex-biased (trend and ratio2 the same sex but detection is both)
    of.write("\n{0} genes with pronounced female-bias:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalFbiasedPronounced,fBiasedXpronounced,expressXsomeWide['X'].sum(),
            fBiasedXpronounced/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            fBiasedApronounced,expressXsomeWide['A'].sum(),fBiasedApronounced/expressXsomeWide['A'].sum()))
    # Output male counts
    # 1) evidence of sex-limited (detection and trend the same sex)
    of.write("\n{0} genes with evidence of male-limited:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalMlimitedEvidence,mLimitedXevidence,expressXsomeWide['X'].sum(),
            mLimitedXevidence/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            mLimitedAevidence,expressXsomeWide['A'].sum(),mLimitedAevidence/expressXsomeWide['A'].sum()))
    # 2) pronounced sex-limited (detection, trend, and ratio2 all the same)
    of.write("\n{0} genes with pronounced male-limited:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalMlimitedPronounced,mLimitedXpronounced,expressXsomeWide['X'].sum(),
            mLimitedXpronounced/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            mLimitedApronounced,expressXsomeWide['A'].sum(),mLimitedApronounced/expressXsomeWide['A'].sum()))
    # 3) evidence of sex-biased (detection is both, and trend is one sex)
    of.write("\n{0} genes with evidence of male-bias:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalMbiasedEvidence,mBiasedXevidence,expressXsomeWide['X'].sum(),
            mBiasedXevidence/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            mBiasedAevidence,expressXsomeWide['A'].sum(),mBiasedAevidence/expressXsomeWide['A'].sum()))
    # 4) pronounced sex-biased (trend and ratio2 the same sex but detection is both)
    of.write("\n{0} genes with pronounced male-bias:\n\t{1} out of {2} expressed genes on the X ({3:.5%})".format(
            totalMbiasedPronounced,mBiasedXpronounced,expressXsomeWide['X'].sum(),
            mBiasedXpronounced/expressXsomeWide['X'].sum()))
    of.write("\n\t{0} out of {1} expressed genes on the autosomes ({2:.5%})".format(
            mBiasedApronounced,expressXsomeWide['A'].sum(),mBiasedApronounced/expressXsomeWide['A'].sum()))
    
    # 5) evidence of mixed sex-bias (detection is both and trend is male_and_female)
    of.write("\n{} genes with evidence for male_and_female ({} on X, {} on autosomes):\n{}\n".format(
            totalMFevidence,mfXevidence,mfAevidence,
            annotDF[annotDF['gene_trend_ttest']=="male_and_female"][['fbgn','symbol','xsome']].to_string(index=False)))
    # 6) pronounced mixed sex-bias (detection is both, trend and ratio2 are male_and_female)
    of.write("\n{} genes with pronounced male_and_female ({} on X, {} on autosomes):\n{}\n".format(
            totalMFpronounced,mfXpronounced,mfApronounced,annotDF[(annotDF['gene_ratio2_ttest']=="male_and_female")
            &(annotDF['gene_trend_ttest']=="male_and_female")][['fbgn','symbol','xsome']].to_string(index=False)))

    # Output the male-limited and female-limited genes
    of.write("\n{} genes with evidence for sex-limited ({} with pronounced):\nfemale-limited({} with evidence, {} pronounced)\n{}\n\nmale-limited({} with evidence, {} pronounced)\n{}\n".format(
            totalFlimitedEvidence+totalMlimitedEvidence,totalFlimitedPronounced+totalMlimitedPronounced,
            totalFlimitedEvidence,totalFlimitedPronounced,annotDF[(annotDF['gene_trend_ttest']=="female")&
                    (annotDF['gene_rna']=="fem")][['fbgn','symbol','xsome','gene_trend_ttest','gene_ratio2_ttest']].to_string(index=False),
            totalMlimitedEvidence,totalMlimitedPronounced,annotDF[(annotDF['gene_trend_ttest']=="male")&
                    (annotDF['gene_rna']=="male")][['fbgn','symbol','xsome','gene_trend_ttest','gene_ratio2_ttest']].to_string(index=False)))
    of.close()
  
    
if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()