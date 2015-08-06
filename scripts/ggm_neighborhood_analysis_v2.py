#!/usr/bin/env python

import os
import logging
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle

def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, filemode='w', level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def readData(fname):
    """ Importing a large DOT file is slow. This function will read a pickle
    file if available. If no pickle, then read DOT and create a pickle for next
    time. """

    pname = os.path.splitext(fname)[0] + ".gpickle"
    try:
        # If there is a pickle, unpickle it
        logging.info("Unpickling file")
        nxGraph = nx.Graph(nx.read_gpickle(pname))
    except:
        logging.info("No Pickled file, will import DOT")
        try:
            # No pickle, try the dot file
            logging.info("Importing dot file")
            nxGraph = nx.Graph(nx.read_dot(fname))

            # Make pickle for next time
            logging.info("Pickle graph for later use.")
            nx.write_gpickle(nxGraph,pname)
        except Exception:
            logging.exception("Please provide a DOT formated file.")
    return(nxGraph)

def getNeighbors(nxGraph, target, geneList, splicingFactors):
    """ Search the primary and secondary neighborhoods. Return neighbors and counts. """
    # Pull primary neighbor from target
    primary = nxGraph.neighbors(target)

    # Pull secondary neighbors
    secondary = list()
    for target2 in primary:
        secondary.extend(nxGraph.neighbors(target2))

    # Calculate the number of primary and secondary neighbors that are in my
    # gene list of interest.
    primSet = set(primary)
    secSet = set([x for x in (primary+secondary) if x != target])  # remove the target gene from the secondary neighborhood.
    sexPrim = primSet & set(geneList)
    sexSec = secSet & set(geneList)
    numPrimary = len(primSet)
    numSecondary = len(secSet)
    numPrimSex = len(sexPrim)
    numSecSex = len(sexSec)
    numSex = len(geneList)

    flag_sex = 0
    if target in geneList:
        flag_sex = 1

    flag_splice = 0
    if target in splicingFactors:
        flag_splice = 1

    return([target, flag_sex, flag_splice, numPrimary, numPrimSex, numSecondary, numSecSex, numSex])

def writeHeader(headerList, handle):
    """ Write a header on to the csv file """
    handle.write(','.join([str(x) for x in headerList]) + "\n")

def writeOutput(myOutList, handle):
    """ Write output CSV """
    handle.write(",".join([str(x) for x in myOutList]) + "\n")

if __name__ == "__main__":

    lname = "/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_neighborhood_analysis.log"
    dotname = "/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.dot"
    oname = "/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_neighborhood_analysis.csv"

    # Turn on Logging if option --log was given
    setLogger(lname,logging.INFO)

    # Import Dot File
    mygraph = readData(dotname)

    # Create gene list by pulling all genes that don't start with 'CG'
    logging.info("Creating gene list")
    geneList = [x for x in mygraph.nodes_iter(data=False) if not x.startswith('CG')]

    # I am wanting to highlight splicing factors in my output graph
    logging.info("Creating splicing factor list")
    splicingFactors = ['vir', 'Rbp1_PD', 'Rbp1_PA', 'B52_PI', 'B52_PK', 'B52_PM', 'B52_PN', 'B52_PO', 'B52_PA', 'B52_PB', 'B52_PC',
            'B52_PD', 'B52_PF', 'sqd_PE', 'sqd_PD', 'sqd_PA', 'sqd_PC', 'sqd_PB', 'Psi', 'mub_PL', 'mub_PK', 'mub_PG', 'mub_PA', 'Rm62',
            'snf', 'mub_PM', 'mub_PI', 'mub_PH', 'mub_PE', 'mub_PD', 'mub_PF', 'mub_PC', 'Spf45_PA', 'ps', 'Spf45_PC', 'Spf45_PB']
    
    # Explore the nieghborhood and make CSV table and subplots
    logging.info("Finding neighbors and writing output")
    with open(oname, 'w') as OUT:
        myHeader = ["FG", "FG_sex", "FG_splice", "num_primary", "num_primary_sex", "num_secondary", "num_secondary_sex", "num_possible_sex"]
        writeHeader(myHeader, OUT)

        # Iterate through all genes in Sex Det
        for node in mygraph.nodes_iter(data=False):
            # Calculate primary and secondary nearest neighbors
            myOut = getNeighbors(mygraph, node, geneList, splicingFactors)
            writeOutput(myOut, OUT)

    logging.info("Script Complete")
