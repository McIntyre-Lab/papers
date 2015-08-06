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

def getNeighbors(nxGraph,target,geneList):
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
    secSet = set(secondary)
    sexPrim = primSet & set(geneList)
    sexSec = secSet & set(geneList)
    numPrimary = len(primSet)
    numSecondary = len(secSet)
    numPrimEnrich = len(sexPrim)
    numSecEnrich = len(sexSec)
    numEnrich = numPrimEnrich + numSecEnrich 
    return(primSet, secSet, sexPrim, sexSec, numPrimary, numSecondary, numPrimEnrich, numSecEnrich, numEnrich)

def getFigName(target, odir):
    """ Create an output name for creating all of the neighborhood subgraphs """
    figName = target + "_ggm_subgraph.png"
    oname = os.path.join(odir, figName)
    return(oname)

def graphSubplot(nxGraph, target, primary, secondary, geneList, spliceList, oname):
    """ Plot primary and secondary neighborhood subgraph """
    subsetList = [target] + list(primary) + list(secondary & set(geneList))
    nxSubGraph = nx.Graph(nx.subgraph(nxGraph, subsetList))

    labelDict = dict()
    colorList = list()
    sizeList = list()
    for gene in nxSubGraph.nodes_iter():
        if gene in spliceList:
            # If gene is a splicing factor color it purple if FG otherwise color royal blue
            if gene == target:
                labelDict[gene] = gene
                sizeList.append(3000)
                colorList.append('#d000ff')
            elif gene in geneList:
                labelDict[gene] = gene
                sizeList.append(3000)
                colorList.append('#007ffd')
        elif gene in geneList:
            # If gene is in Sex Det color it red if FG otherwise color light blue
            if gene == target:
                labelDict[gene] = gene
                sizeList.append(3000)
                colorList.append('red')
            else:
                labelDict[gene] = gene
                sizeList.append(3000)
                colorList.append('#aeeeee')
        else:
            # If gene is not in Sex Det make it small and color it green if 'CG7705' otherwise color grey 
            if gene == "CG7705":
                labelDict[gene] = ''
                sizeList.append(300)
                colorList.append('#00f554')
            elif gene == "CG1406":
                labelDict[gene] = ''
                sizeList.append(300)
                colorList.append('#FAFF00')
            else:
                labelDict[gene] = ''
                sizeList.append(300)
                colorList.append('grey')

    # Generate layout using force spring algorithm
    pos = nx.spring_layout(nxSubGraph, iterations=200, k=.2)

    # Plot Subgraph
    fig = plt.figure(figsize=(20,20),dpi=300)
    nx.draw_networkx_nodes(nxSubGraph,pos,node_size=sizeList,node_color=colorList)
    nx.draw_networkx_labels(nxSubGraph,pos,labelDict)
    nx.draw_networkx_edges(nxSubGraph,pos)
    plt.axis('off')
    fig.savefig(oname, format='png', bbox_inches='tight')
    plt.close()

def writeHeader(handle):
    """ Write a header on to the csv file """
    handle.write("gene,primary,secondary,sexPrimary, sexSecondary, num_primary, num_secondary, num_primaryEnrich, num_secondaryEnrich, num_enrich\n")

def writeOutput(gene, primary, secondary, sexPrim, sexSec, numPrimary, numSecondary, numPrimEnrich, numSecEnrich, numEnrich, handle):
    """ Write output CSV """
    prim = '|'.join([str(x) for x in primary])
    sec = '|'.join([str(x) for x in secondary])
    pSet = '|'.join([str(x) for x in sexPrim])
    sSet = '|'.join([str(x) for x in sexSec])
    handle.write("{0},{1},{2},{3},{4},{5},{6},{7},{8},{9}\n".format(gene, prim, sec, pSet, sSet, numPrimary, numSecondary, numPrimEnrich, numSecEnrich, numEnrich))
    

if __name__ == "__main__":

    lname = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2_neighbor_table.log'
    dname = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2.dot'
    oname = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/dsrp_ggm_isoforms_FDR2_neighbor_table.csv'
    odir = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/ggm/subgraphs'

    # Turn on Logging if option --log was given
    setLogger(lname,logging.INFO)

    # Import Dot File
    mygraph = readData(dname)

    # Create gene list by pulling all genes that don't start with 'CG'
    logging.info("Creating gene list")
    geneList = [x for x in mygraph.nodes_iter(data=False) if not x.startswith('CG')]

    # I am wanting to highlight splicing factors in my output graph
    logging.info("Creating splicing factor list")
    splicingFactors = ['vir', 'Rbp1_PD', 'Rbp1_PA', 'B52_PI', 'B52_PK',
            'B52_PM', 'B52_PN', 'B52_PO', 'B52_PA', 'B52_PB', 'B52_PC',
            'B52_PD', 'B52_PF', 'sqd_PE', 'sqd_PD', 'sqd_PA', 'sqd_PC',
            'sqd_PB', 'Psi', 'mub_PL', 'mub_PK', 'mub_PG', 'mub_PA', 'Rm62',
            'snf', 'mub_PM', 'mub_PI', 'mub_PH', 'mub_PE', 'mub_PD', 'mub_PF',
            'mub_PC', 'Spf45_PA', 'ps', 'Spf45_PC', 'Spf45_PB']
    
    # Explore the nieghborhood and make CSV table and subplots
    logging.info("Finding neighbors and writing output")
    with open(oname, 'w') as OUT:
        writeHeader(OUT)

        # Iterate through all genes in Sex Det
        for gene in geneList:
            fName = getFigName(gene, odir)

            # Calculate primary and secondary nearest neighbors
            primary, secondary, sexPrim, sexSec, numPrimary, numSecondary, numPrimEnrich, numSecEnrich, numEnrich = getNeighbors(mygraph, gene, geneList)
            writeOutput(gene, primary, secondary, sexPrim, sexSec, numPrimary, numSecondary, numPrimEnrich, numSecEnrich, numEnrich, OUT)

            # Create neighborhood subplots so neighborhoods can be visuzlized
            graphSubplot(mygraph, gene, primary, secondary, geneList, splicingFactors, fName)

    logging.info("Script Complete")
