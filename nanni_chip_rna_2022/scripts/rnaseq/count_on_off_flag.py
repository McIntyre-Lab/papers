#!/usr/bin/env python

import argparse
import pandas as pd

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Count on/off flags")

    # Input data
    parser.add_argument("-f", "--flags", dest="inFile", required=True, help="Input on/off flag CSV file")

    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output text file name of counts")

    args = parser.parse_args()
    return args

def main():
    # Get input flag file
    flagDF = pd.read_csv(args.inFile)
    
    # Open output file
    of = open(args.outFile,'w')
    
    # Get crosstabs of male/female flags
    of.write("\n"+str(pd.crosstab(flagDF.filter(regex="_f_on").iloc[:,0], 
                             flagDF.filter(regex="_m_on").iloc[:,0]))+"\n")
    of.write("\n"+str(pd.crosstab(flagDF.filter(regex="_f_allsID_on").iloc[:,0], 
                             flagDF.filter(regex="_m_allsID_on").iloc[:,0]))+"\n")
    of.write("\n"+str(pd.crosstab(flagDF.filter(regex="_f_allsID_on").iloc[:,0], 
                             flagDF.filter(regex="_m_allsID_off").iloc[:,0]))+"\n")
    of.write("\n"+str(pd.crosstab(flagDF.filter(regex="_f_allsID_off").iloc[:,0], 
                             flagDF.filter(regex="_m_allsID_on").iloc[:,0]))+"\n")
    of.close()

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

