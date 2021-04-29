#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Count DE Fold Change Consistency")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Input CSV of combiend DE tappAS results (with or without GO)")

    # Output data
    parser.add_argument("-g", "--out-group", dest="outGroup", required=True, help="Output CSV file of group FC direction gene counts")
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV of tappAS results with new flags")

    args = parser.parse_args()
    return args

def main():
    # Get input file
    inDF = pd.read_csv(args.inFile)
    
    # Make extra flags for GO
    inDF['flag_detect_DE_NC338_only'] = np.where((inDF['flag_detect_DE_B73']==0)&
                                                 (inDF['flag_detect_DE_C123']==0)&
                                                 (inDF['flag_detect_DE_Mo17']==0)&
                                                 (inDF['flag_detect_DE_Hp301']==0)&
                                                 (inDF['flag_detect_DE_NC338']==1),1,0)
    inDF['flag_detect_DE_all_noB73'] = np.where((inDF['flag_detect_DE_B73']==0)&
                                                 (inDF['flag_detect_DE_C123']==1)&
                                                 (inDF['flag_detect_DE_Mo17']==1)&
                                                 (inDF['flag_detect_DE_Hp301']==1)&
                                                 (inDF['flag_detect_DE_NC338']==1),1,0)
    
    # Get groups of FC direction
    for genotype in ["B73","C123","Hp301","Mo17","NC338"]:
        inDF[genotype+'_FC_direction'] = np.where(inDF['Log2FC_'+genotype]>0,"UP",np.where(
            inDF['Log2FC_'+genotype]<0,"DOWN","zero"))
    
    # Count number of genes DE in at least one genotype
    deDF = inDF[(inDF['flag_detect_DE_B73']==1)|(inDF['flag_detect_DE_C123']==1)|
                (inDF['flag_detect_DE_Mo17']==1)|(inDF['flag_detect_DE_Hp301']==1)|
                (inDF['flag_detect_DE_NC338']==1)].copy()
    print("{} gene DE in at least one genotype".format(len(deDF)))
    
    # Count number of genes with consistent up FC
    cols = [c for c in deDF if "Log2FC" in c]
    upDF = deDF[(deDF[cols] > 0).all(axis=1)]
    downDF = deDF[(deDF[cols] < 0).all(axis=1)]
    print("{0:.2%} ({1} genes total, {2} up, {3} down) have consistent direction of FC".format(
        (len(upDF)+len(downDF))/len(deDF),len(upDF)+len(downDF),len(upDF),len(downDF)))
    
    # Count groups of FC direction
    FCcols = [c for c in deDF.columns if "FC_direction" in c]
    groupDF = deDF.groupby(FCcols)['gene_id'].count().reset_index().sort_values(['gene_id'],ascending=[False]).rename(columns={'gene_id':'num_gene'})
    
    groupDF.to_csv(args.outGroup,index=False)
    
    # Make flags for groups of interest
    inDF['flag_B73Up_restDown'] = np.where(((inDF['flag_detect_DE_B73']==1)|
                                            (inDF['flag_detect_DE_C123']==1)|
                                            (inDF['flag_detect_DE_Mo17']==1)|
                                            (inDF['flag_detect_DE_Hp301']==1)|
                                            (inDF['flag_detect_DE_NC338']==1))&
                                           (inDF['B73_FC_direction']=="UP")&
                                           (inDF['C123_FC_direction']=="DOWN")&
                                           (inDF['Hp301_FC_direction']=="DOWN")&
                                           (inDF['Mo17_FC_direction']=="DOWN")&
                                           (inDF['NC338_FC_direction']=="DOWN"),1,0)
    inDF['flag_B73Zero_restUp'] = np.where(((inDF['flag_detect_DE_B73']==1)|
                                            (inDF['flag_detect_DE_C123']==1)|
                                            (inDF['flag_detect_DE_Mo17']==1)|
                                            (inDF['flag_detect_DE_Hp301']==1)|
                                            (inDF['flag_detect_DE_NC338']==1))&
                                           (inDF['B73_FC_direction']=="zero")&
                                           (inDF['C123_FC_direction']=="UP")&
                                           (inDF['Hp301_FC_direction']=="UP")&
                                           (inDF['Mo17_FC_direction']=="UP")&
                                           (inDF['NC338_FC_direction']=="UP"),1,0)
    inDF['flag_B73Down_restUp'] = np.where(((inDF['flag_detect_DE_B73']==1)|
                                            (inDF['flag_detect_DE_C123']==1)|
                                            (inDF['flag_detect_DE_Mo17']==1)|
                                            (inDF['flag_detect_DE_Hp301']==1)|
                                            (inDF['flag_detect_DE_NC338']==1))&
                                           (inDF['B73_FC_direction']=="DOWN")&
                                           (inDF['C123_FC_direction']=="UP")&
                                           (inDF['Hp301_FC_direction']=="UP")&
                                           (inDF['Mo17_FC_direction']=="UP")&
                                           (inDF['NC338_FC_direction']=="UP"),1,0)
    inDF['flag_NC338Down_restUp'] = np.where(((inDF['flag_detect_DE_B73']==1)|
                                            (inDF['flag_detect_DE_C123']==1)|
                                            (inDF['flag_detect_DE_Mo17']==1)|
                                            (inDF['flag_detect_DE_Hp301']==1)|
                                            (inDF['flag_detect_DE_NC338']==1))&
                                           (inDF['B73_FC_direction']=="UP")&
                                           (inDF['C123_FC_direction']=="UP")&
                                           (inDF['Hp301_FC_direction']=="UP")&
                                           (inDF['Mo17_FC_direction']=="UP")&
                                           (inDF['NC338_FC_direction']=="DOWN"),1,0)
    # Output new flags
    inDF.to_csv(args.outFile,index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
