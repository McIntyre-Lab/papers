#!/usr/bin/env python
 
import argparse
import pandas as pd
import numpy as np
import os

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Flag monoexon PacBio isoforms and reference transcripts in and output counts")

    # Input arguments
    parser.add_argument("-c", "--classification", dest="inFile", required=True, help="PacBio classification file from SQANTI QC")

    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output file for frequencies")

    args = parser.parse_args()
    return args


def main():
    # Check that input files exits
    if not os.path.isfile(args.inFile) :
        raise FileNotFoundError

    # Get input classification file
    classDF = pd.read_csv(args.inFile, sep="\t", low_memory=False)
    
    # Make flags based on exon numbers
    classDF['flag_PB_monoexon'] = np.where(classDF['exons']==1,1,0)
    classDF['flag_ref_monoexon'] = np.where(classDF['ref_exons']==1,1,0)
    
    # Get frequencies of flags
    f = open(args.outFile,'w')
    f.write(str(pd.crosstab(classDF['flag_PB_monoexon'],classDF['flag_ref_monoexon']))+"\n")
    f.close()
    
    # Output flag classification file
    prefix = args.inFile.strip(".txt").strip(".tsv")
    classDF.to_csv(prefix+"_flag_mono.txt",sep='\t',index=False)
    
    # Looking at the isoforms that are monoexon matches to multiexon reference
    mono2multi = classDF[(classDF['flag_PB_monoexon']==1) & (classDF['flag_ref_monoexon']==0)].copy()
    print("\nMonoexon PB matched to multiexon reference:\n")
    print(mono2multi.groupby('structural_category')['isoform'].count().to_string()+"\n")
    print(mono2multi.groupby('subcategory')['isoform'].count().to_string()+"\n")
    print("PB length\n"+mono2multi['length'].describe().to_string()+"\n")
    print("Percent A downstream TTS\n"+mono2multi['perc_A_downstream_TTS'].describe().to_string()+"\n")
    mono2multi['PB2ref_length_ratio'] = mono2multi['length']/mono2multi['ref_length']
    print("PB2ref length ratio\n"+mono2multi['PB2ref_length_ratio'].describe().to_string()+"\n")
    
    # Looking at the isoforms that are monoexon matches to monoexon reference
    mono2mono = classDF[(classDF['flag_PB_monoexon']==1) & (classDF['flag_ref_monoexon']==1)].copy()
    print("\nMonoexon PB matched to monoexon reference:\n")
    print(mono2mono.groupby('structural_category')['isoform'].count().to_string()+"\n")
    print(mono2mono.groupby('subcategory')['isoform'].count().to_string()+"\n")
    print("PB length\n"+mono2mono['length'].describe().to_string()+"\n")
    print("Percent A downstream TTS\n"+mono2mono['perc_A_downstream_TTS'].describe().to_string()+"\n")
    mono2mono['PB2ref_length_ratio'] = mono2mono['length']/mono2mono['ref_length']
    print("PB2ref length ratio\n"+mono2mono['PB2ref_length_ratio'].describe().to_string()+"\n")
    
    # Looking at the isoforms that are multiexon matches to multiexon reference
    multi2multi = classDF[(classDF['flag_PB_monoexon']==0) & (classDF['flag_ref_monoexon']==0)].copy()
    print("\nMultiexon PB matched to multiexon reference:\n")
    print(multi2multi.groupby('structural_category')['isoform'].count().to_string()+"\n")
    print(multi2multi.groupby('subcategory')['isoform'].count().to_string()+"\n")
    print("PB length\n"+multi2multi['length'].describe().to_string()+"\n")
    print("Percent A downstream TTS\n"+multi2multi['perc_A_downstream_TTS'].describe().to_string()+"\n")
    multi2multi['PB2ref_length_ratio'] = multi2multi['length']/multi2multi['ref_length']
    print("PB2ref length ratio\n"+multi2multi['PB2ref_length_ratio'].describe().to_string()+"\n")
     
if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()

