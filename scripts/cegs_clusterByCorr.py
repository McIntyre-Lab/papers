#!/usr/bin/env python
import argparse 
import os.path
import logging
import numpy as np
import scipy as sc
import pandas as pd

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Takes a dataset and collapses columns based on correlation coefficients.")
    parser.add_argument("-i", dest="fname", action='store', required=True, help="Input CSV file [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Name of the output CSV file [Required]")
    parser.add_argument("-t", dest="tname", action='store', required=True, help="Name of the output table containing the list of collapse isoform ids [Required]")
    parser.add_argument("-c", dest="cutoff", action='store', required=False, default=.8, help="Cut-off value to group by [Default 0.8]")
    parser.add_argument("-g", "--log", dest="log", action='store', required=False, help="Path and name of log file [Optional]") 
    args = parser.parse_args()
    #args = parser.parse_args(['-i', '/home/jfear/tmp/b52_test.csv', '-o', '/home/jfear/tmp/b52_merged.csv', '-t', '/home/jfear/tmp/b52_merged.tbl', '-g', '/home/jfear/tmp/b52_merged.log'])
    return(args)

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, level=loglevel, filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')

def importDataFile(fname):
    """ Import the data file in a CSV format. Note assuming that the sample
    identifiers are present 'line'. Once file is in, figure out what gene this
    is. 
    """
    logging.info("Importing file: {0}".format(fname))

    # Read in file to pandas dataframe and figure how many rows and columns there are.
    df = pd.read_csv(fname, index_col=['line'])
    rows, cols = df.shape

    # Using the headers grab gene information
    name = df.columns.values.tolist()[0].rsplit('_',1)[0]
    logging.info("Current gene symbol: {0}".format(name))
    logging.info("{0} columns and {1} rows were successfully imported".format(cols, rows))
    return(df, name)

def getMaskList(cutoff, myCorr):
    """ Create a boolean mask of the correlation matrix, where correlation
    exceeds the given cutoff (.8) by default. For each row of the mask, create
    a list of coordinates in the correlation matrix that are True set. 
    """
    logging.info("Identifying regions that pass the cutoff filter.")
    mask = myCorr >= float(cutoff)
    myList = []
    for index, row in mask.iterrows():
        myList.append(np.where(row)[0])

    # Add warning to log if correlations are pretty high, but not quite the
    # cutoff. So if the value is below our cutoff, but greater than 0.75 then I
    # will flag with a warning.
    mask2 = np.logical_not(mask)
    mask2 &= myCorr > .75
    if mask2.any().any():
        logging.warn("There are some large correlation values that did not reach the cutoff. You may want to hand curate these: \n{0}".format(mask2.to_string()))
    return(myList)

def mergeReduce(myList):
    """ Takes a list of masked coordinates, converts these to a list of sets of
    coordinates. Then keeps iterating over the list combining overlapping sets.
    Returns list of combined sets.  
    """
    logging.info("Reducing overlapping isoforms.")
    sets = [set(lst.tolist()) for lst in myList]
    merged = 1
    while merged:
        merged = 0
        results = []
        while sets:
            common, rest = sets[0], sets[1:]
            sets = []
            for x in rest:
                if x.isdisjoint(common):
                    sets.append(x)
                else:
                    merged = 1
                    # if they match take the union
                    common |= x
            results.append(common)
        sets = results
    return(sets)

def alphaList():
    import string
    return(list(string.ascii_uppercase))

def createCollapseTable(df, currSet):
    name = [i.rsplit('_', 1)[1] for i in df.columns.values.tolist()]
    myGroup = ';'.join([name[i] for i in currSet])
    return(myGroup)

def createNewIso(df, geneName, mySet, myAlpha):
    """ Average together isoforms that were identified as being highly
    correlated. Then return a dictionary of results.
    """
    logging.info("Averaging overlapping isoforms.")
    myDict = dict()
    myTable = dict()
    for ind in xrange(len(mySet)):
        if len(mySet) == 1:
            name = geneName
            myTable[name] = name
        else:
            name = geneName + '_' + myAlpha[ind]
            myTable[name] = createCollapseTable(df, mySet[ind])
        value = df.iloc[:, mySet[ind]].apply(np.mean, axis=1)
        myDict[name] = value.round(decimals=10)
    return(myDict,myTable)

def main(args):
    # Get file name and import into a data frame
    df, gname = importDataFile(args.fname)

    # Calculate the correlation and figure out where all of the correlation coefficients is > 0.8
    logging.info("Calculating correlation structure.")
    myCorr = df.corr()
    logging.info("Correlation matrix: \n{0}".format(myCorr.to_string()))

    # Get a list of locations that the mask was true
    myList = getMaskList(args.cutoff, myCorr)

    # Merge the list of mask locations to get sets of isoforms grouped together
    mySet = mergeReduce(myList)

    # Now create new "isoform" groups with an alpha numeric index append it to the gene name
    myAlpha = alphaList()
    newIso, isoTable = createNewIso(df, gname, mySet, myAlpha)

    # Combine these into a data frame and write out
    logging.info("Exporting new dataset to: {0}".format(args.oname))
    newDf = pd.DataFrame(newIso)
    newDf.to_csv(args.oname)

    # Create Table with collapsed isoform names 
    logging.info("Exporting collapsed isoform names: {0}".format(args.tname))
    with open(args.tname, 'w') as OUT:
        OUT.write("name,isoform\n")
        for key, value in isoTable.iteritems():
            OUT.write("{0},{1}\n".format(key, value))

if __name__ == '__main__':
    # Turn on Logging if option -g was given
    args = getOptions()
    if args.log:
        setLogger(args.log,logging.INFO)
    else:
        setLogger(os.devnull,logging.INFO)

    # Run Main part of the script
    logging.info("Starting script.")
    main(args)
    logging.info("Script complete.")

