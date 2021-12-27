#!/usr/bin/env python3

import pandas as pd
import argparse

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Takes in sex specifidc flag file (*_flag.csv from count_detected_features_08avn.py) and calculates kappa agreement")

    # Input data
    parser.add_argument("-i", "--input-flags", dest="inFlags", required=True, help="Input CSV of flags, has featureID and a column of each detection (flag_sex_antibody_on)")

    # Output data
    parser.add_argument("-o", "--output", dest="outPrefix", required=True, help="Prefix for output files")

    args = parser.parse_args()
    return args

def main():

    # Get input flag file
    flagDF = pd.read_csv(args.inFlags)

    # Output flag crosstabs and kappa agreement values
    f = open(args.outPrefix+"_flag_crosstabs.csv",'w')
    k = open(args.outPrefix+"_flag_kappas.csv",'w')
    k.write("featureType,K4_po,K27_po,f_po,m_po,K4_pe,K27_pe,f_pe,m_pe,K4_kappa,K27_kappa,f_kappa,m_kappa")
    for FEAT in flagDF['featureType'].unique():
        tempDF = flagDF[flagDF['featureType']==FEAT]
        tempOpen = pd.crosstab(tempDF['flag_f_K4_on'],tempDF['flag_m_K4_on'])
        tempClose = pd.crosstab(tempDF['flag_f_K27_on'],tempDF['flag_m_K27_on'])
        tempFemale = pd.crosstab(tempDF['flag_f_K4_on'],tempDF['flag_f_K27_on'])
        tempMale = pd.crosstab(tempDF['flag_m_K4_on'],tempDF['flag_m_K27_on'])
        f.write("\n"+str(FEAT)+"\n"+tempOpen.to_string()+"\n")
        f.write("\n"+tempClose.to_string()+"\n\n")
        K4po = (tempOpen.iloc[0,0] + tempOpen.iloc[1,1])/ tempOpen.values.sum()
        K27po = (tempClose.iloc[0,0] + tempClose.iloc[1,1])/ tempClose.values.sum()
        Fpo = (tempFemale.iloc[0,0] + tempFemale.iloc[1,1])/ tempFemale.values.sum()
        Mpo = (tempMale.iloc[0,0] + tempMale.iloc[1,1])/ tempMale.values.sum()
        K4pe = ((tempOpen.iloc[[0]].values.sum()/tempOpen.values.sum())*(tempOpen[0].values.sum()/tempOpen.values.sum()))+((tempOpen.iloc[[1]].values.sum()/tempOpen.values.sum())*(tempOpen[1].values.sum()/tempOpen.values.sum()))
        K27pe = ((tempClose.iloc[[0]].values.sum()/tempClose.values.sum())*(tempClose[0].values.sum()/tempClose.values.sum()))+((tempClose.iloc[[1]].values.sum()/tempClose.values.sum())*(tempClose[1].values.sum()/tempClose.values.sum()))
        Fpe = ((tempFemale.iloc[[0]].values.sum()/tempFemale.values.sum())*(tempFemale[0].values.sum()/tempFemale.values.sum()))+((tempFemale.iloc[[1]].values.sum()/tempFemale.values.sum())*(tempFemale[1].values.sum()/tempFemale.values.sum()))
        Mpe = ((tempMale.iloc[[0]].values.sum()/tempMale.values.sum())*(tempMale[0].values.sum()/tempMale.values.sum()))+((tempMale.iloc[[1]].values.sum()/tempMale.values.sum())*(tempMale[1].values.sum()/tempMale.values.sum()))
        K4kappa = (K4po - K4pe)/(1 - K4pe)
        K27kappa = (K27po - K27pe)/(1 - K27pe)
        Fkappa = (Fpo - Fpe)/(1 - Fpe)
        Mkappa = (Mpo - Mpe)/(1 - Mpe)
        k.write("\n{},{},{},{},{},{},{},{},{},{},{},{},{}".format(FEAT,K4po,K27po,Fpo,Mpo,K4pe,K27pe,Fpe,Mpe,K4kappa,K27kappa,Fkappa,Mkappa))
    f.close()
    k.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
