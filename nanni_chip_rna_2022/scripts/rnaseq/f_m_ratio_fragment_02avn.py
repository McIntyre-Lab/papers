#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Get ratios of female/male APN")

    # Input data
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Input coverage counts CSV file")

    # Output data
    parser.add_argument("-o", "--outPrefix", dest="outPrefix", required=True, help="Output prefix")

    args = parser.parse_args()
    return args

def main():
    # Get input coverage file
    covDF = pd.read_csv(args.inFile)
    
    # Get average of male and female uq normalized samples
    # Where avg is 0 set to 10^-6
    fcols = [c for c in covDF.columns if ("_f_" in c) and ("_apn_uq_ff" in c)]
    mcols = [c for c in covDF.columns if ("_m_" in c) and ("_apn_uq_ff" in c)]
    covDF['avg_f_apn_uq_ff'] = covDF[fcols].mean(axis=1).replace(0,10**-6)
    covDF['avg_m_apn_uq_ff'] = covDF[mcols].mean(axis=1).replace(0,10**-6)
    
    # Get ratio of the average female / average male
    covDF['ratio_avg_f_m_apn_uq_ff'] = covDF['avg_f_apn_uq_ff']/covDF['avg_m_apn_uq_ff']
    
    # Print ratio value counts in ranges to log
    f = open(args.outPrefix+"_counts.txt",'w')
    bins=[covDF['ratio_avg_f_m_apn_uq_ff'].min(),0.2,0.25,0.33,0.5,1,2,3,4,5,covDF['ratio_avg_f_m_apn_uq_ff'].max()]
    f.write("ratio\t\tcount\n")
    f.write(covDF['ratio_avg_f_m_apn_uq_ff'].value_counts(bins=bins,sort=False).to_string()+"\n")

    # Flag male (=<0.5) and female (>=2)
    covDF['flag_f_bias_RNA_apn_uq_ff'] = np.where(covDF['ratio_avg_f_m_apn_uq_ff']>=2,1,0)
    covDF['flag_m_bias_RNA_apn_uq_ff'] = np.where(covDF['ratio_avg_f_m_apn_uq_ff']<=0.5,1,0)
    covDF['flag_no_bias_RNA_apn_uq_ff'] = np.where((covDF['ratio_avg_f_m_apn_uq_ff']>0.5) & (covDF['ratio_avg_f_m_apn_uq_ff']<2),1,0)
    f.write("\nsum(flag_f_bias_RNA_apn_uq_ff) = "+str(covDF['flag_f_bias_RNA_apn_uq_ff'].sum())+"\n")
    f.write("sum(flag_m_bias_RNA_apn_uq_ff) = "+str(covDF['flag_m_bias_RNA_apn_uq_ff'].sum())+"\n")
    f.write("sum(flag_no_bias_RNA_apn_uq_ff) = "+str(covDF['flag_no_bias_RNA_apn_uq_ff'].sum())+"\n")
    f.write("total = "+str(len(covDF))+"\n")
    f.close()

    # Output flags
    covDF[['featureID','avg_f_apn_uq_ff','avg_m_apn_uq_ff','ratio_avg_f_m_apn_uq_ff','flag_f_bias_RNA_apn_uq_ff','flag_m_bias_RNA_apn_uq_ff',
              'flag_no_bias_RNA_apn_uq_ff']].to_csv(args.outPrefix+"_sex_bias.csv", index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

