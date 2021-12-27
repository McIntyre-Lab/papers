
#!/usr/bin/env python
 

import argparse
import pandas as pd
import numpy as np
import sqlite3
import os
import sys

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Take in ChIP-seq feature Detection Above Input (DAI) flag file and create BED file of average corrected APN values (0 or features not detected)")

    # Input arguments
    parser.add_argument("-i", "--input", dest="inFile", required=True, help="Feature DAI flag CSV file")
    parser.add_argument("--fragment", dest="inFrag", required=True, help="Fragment BED file corresponding to the fragments analyzed in the flag file")
    parser.add_argument("--fusion", dest="inFus", required=True, help="Fusion BED file corresponding to the fusions analyzed in the flag file")
    parser.add_argument("--intron", dest="inIntr", required=True, help="Intron BED file corresponding to the introns analyzed in the flag file")
    parser.add_argument("-d", "--db-file", dest="tempDB", required=False, help="Temp SQL database file to reduce memory required", default=":memory:")

    # Output arguments
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output BED file name")

    args = parser.parse_args()
    return args


def main():
    # Check that input files exits
    if not os.path.isfile(args.inFile) or not os.path.isfile(args.inFrag) or not os.path.isfile(args.inFus) or not os.path.isfile(args.inIntr) :
        raise FileNotFoundError

    # Get input flag file
    flagDF = pd.read_csv(args.inFile)
    
    # Get feature files
    fragDF = pd.read_csv(args.inFrag, sep="\t", names=['chr','start','end','featureID'])
    fusDF = pd.read_csv(args.inFus, sep="\t", names=['chr','start','end','featureID'])
    intronDF = pd.read_csv(args.inIntr, sep="\t", names=['chr','start','end','featureID'])
    
    # Stack feature files
    eaDF = pd.concat([fragDF,fusDF,intronDF], ignore_index=True)
    
    # Connect to SQL database
    con = sqlite3.connect(args.tempDB)
    cur = con.cursor()
    
    # Merge feature coordinates with flags
    flagDF.to_sql("flags", con, if_exists="replace")
    eaDF.to_sql("features", con, if_exists="replace")
    cur.execute("CREATE TABLE merge AS SELECT in1.*, in2.chr, in2.start, in2.end "
                "FROM flags in1 LEFT JOIN features in2 "
                "ON in1.featureID = in2.featureID ;")
    mergeDF = pd.read_sql("SELECT * FROM merge", con).drop(columns=['index'])
    
    # Parse TSS, 3UTR, and 5UTR names to get coordinates
    mergeDF['chr'] = np.where((mergeDF['featureID'].str.contains("intergenic")) | 
           (mergeDF['featureID'].str.contains("3UTR")) |
           (mergeDF['featureID'].str.contains("5UTR")) |
           (mergeDF['featureID'].str.contains("TSS")), 
           mergeDF['featureID'].str.split("_").str[1:-2].str.join("_"), mergeDF['chr'])
    mergeDF['start'] = np.where((mergeDF['featureID'].str.contains("intergenic")) | 
           (mergeDF['featureID'].str.contains("3UTR")) |
           (mergeDF['featureID'].str.contains("5UTR")) |
           (mergeDF['featureID'].str.contains("TSS")), 
           mergeDF['featureID'].str.split("_").str[-2], mergeDF['start'])
    mergeDF['end'] = np.where((mergeDF['featureID'].str.contains("intergenic")) | 
           (mergeDF['featureID'].str.contains("3UTR")) |
           (mergeDF['featureID'].str.contains("5UTR")) |
           (mergeDF['featureID'].str.contains("TSS")), 
           mergeDF['featureID'].str.split("_").str[-1], mergeDF['end'])
    
    # Check for proper parsing (no NA values left in dataframe)
    if len(mergeDF[mergeDF.isnull().any(axis=1)]) > 0:
        print("WARNING : Incorrect parsing of feature coordinates - check input featureID format")
        sys.exit(1)
    
    # Set 0 apnFF values to NaN (so they will not be included in averages)
    mergeDF[['rep1_apnFF','rep2_apnFF','rep3_apnFF']] = mergeDF[['rep1_apnFF','rep2_apnFF','rep3_apnFF']].replace(0,np.nan)
    
    # Average apn values - where peak not detected (presence = 0) set average = 0
    mergeDF['avg_apnFF'] = np.where(mergeDF['flag_presence']==0,0,
           mergeDF[['rep1_apnFF','rep2_apnFF','rep3_apnFF']].mean(axis=1))
    
    # Print BED output file
    mergeDF[['chr','start','end','avg_apnFF','featureID']].to_csv(args.outFile,sep="\t",index=False,header=False)


if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()
