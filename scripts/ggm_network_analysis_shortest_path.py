#!/usr/bin/env python

import os
import logging
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pickle
import itertools

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

def writeHeader(handle):
    """ Write a header on to the csv file """
    handle.write("gene1,gene2,shortest_distance\n")

def writeOutput(gene1, gene2, dist, handle):
    handle.write("{0},{1},{2}\n".format(gene1, gene2, dist))
    
if __name__ == "__main__":

    dname = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.dot'
    oname = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2_shortest_path_table.csv'
    lname = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2_shortest_path_table.log'

    # Turn on Logging if option --log was given
    setLogger(lname,logging.INFO)

    # Import Dot File
    mygraph = readData(dname)

    # Create gene list by pulling all genes that don't start with 'CG'
    logging.info("Creating gene list")
    geneList = [x for x in mygraph.nodes_iter(data=False) if not x.startswith('CG')]

    # Iterate through all of my genes of interest output their primary and
    # secondary neighbors
    logging.info("Finding neighbors and writing output")
    with open(oname, 'w') as OUT:
        writeHeader(OUT)
        for permut in itertools.combinations(geneList,2):
            dist = nx.shortest_path_length(mygraph, permut[0], permut[1])
            writeOutput(permut[0],permut[1],dist,OUT)

    logging.info("Script Complete")
