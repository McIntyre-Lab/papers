#!/usr/bin/env python

import argparse
import logging
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="This script takes SAS SEM output and graphs it a given network structure.")
    parser.add_argument("--pos", dest="pname", action='store', required=True, help="A CSV file containing 'NodeID,x-pos,y-pos' [Required].")
    parser.add_argument("--sem", dest="sname", action='store', required=True, help="A CSV file containing 'NodeID,x-pos,y-pos' [Required].")
    parser.add_argument("--out", dest="out", action='store', required=True, help="Output Png File", metavar="OUT_FILE")
    parser.add_argument("--log", dest="log", action='store', required=False, help="Log File", metavar="LOG_FILE") 
    args = parser.parse_args()
    #args = parser.parse_args(['--pos', '/home/jfear/mclab/cegs_sem_sd_paper/documentation/networks/sex_det_node_positions_w_isoforms.csv', '--sem', '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/sem/unconstrained_estimates.csv', '--out', '/home/jfear/tmp/bob.png', '--log', '/home/jfear/tmp/bob.log'])
    return(args)


def setLogger(fname,loglevel):
    """ Function to handle error logging """
    logging.basicConfig(filename=fname, filemode='w', level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')

def getNodePos(pname):
    """ Import node positions from CSV file and put information into a
    dictionary where key is the node and the value is a numpy array with x,
    y-coordinates """
    logging.info("Reading Positions File")
    mypos = dict()
    with open(pname, 'r') as CSV:
        for row in CSV:
            cols = row.rstrip().split(',')
            mypos[cols[0]] = np.array([float(cols[1]),float(cols[2])])
    return(mypos)

def getSemResults(sname):
    """ Import SEM results from CSV.  Dataset is created by using ODS to output
    PATHListStd and PATHCovVarsStd.  The 'arrow' variable was droped from
    PATHListStd, and the two datasets were set together. 
    
    Headers are: SpecType,Var1,Var2,Parameter,Estimate,StdErr,tValue

    SpecType is 9 for path from Var1 -> Var2
    SpecType is -14 for covaraince between Var1 <-> Var2

    mysem is a list of significant edges (critValue >= 1.96) in the form of tuples containing (Var1, Var2, edgeAttributes)
    expNodes is a list of nodes that were expressed and modeled.
    """
    logging.info("Reading SEM output File")
    with open(sname, 'r') as CSV:
        CSV.next()  # skip header
        mysem = list()
        nodes = list()
        for row in CSV:
            cols = row.rstrip().split(',')
            SpecType, Var1, Var2, Parameter, Estimate, StdErr, tValue = cols
            nodes.extend((Var1,Var2))

            try:
                # Only add edge if t-value is significant (alpha = 0.05, critVal = 1.96) 
                if abs(float(tValue)) >= 1.96:
                    # If an estimate is negative, then make the edge red to show
                    # repression.
                    if float(Estimate) < 0:
                        edgeColor = 'red'
                    else:
                        edgeColor = 'black'

                    if SpecType == '9':
                        edgeDir=''
                    elif SpecType == '-14':
                        edgeDir='both'
                    else:
                        logging.error("SpecType has been hard coded to be 9 or -14 in this script, but it was %s. Please check SAS output or change this script." % SpecType)
                        raise TypeError

                    if edgeDir:
                        # To get bi-directional arrows for covarainces, I am added to edges in opposite directions
                        edgeAttribute = {'color': edgeColor, 'dir': edgeDir, 'label': Estimate}
                        mytup = (Var1, Var2, edgeAttribute)
                        mysem.append(mytup)
                        mytup = (Var2, Var1, edgeAttribute)
                        mysem.append(mytup)
                    else:
                        edgeAttribute = {'color': edgeColor, 'label': Estimate}
                        mytup = (Var1, Var2, edgeAttribute)
                        mysem.append(mytup)
            except:
                # Constrained edges will have a '_' for estimate and tValue.
                # This will cause an exception in the above if statement. Catch
                # the exception and add a note to the log.
                logging.warn("{0} and {1} have a constrained Edge in SEM output".format(Var1,Var2))

    expNodes = set(nodes)
    return(mysem,expNodes)

def createGraph(pos, sem):
    """ Create a networkx graph object """
    logging.info("Building Network Object")
    G = nx.MultiDiGraph() # Using multi so I can graph covariances as double headed arrorws
    G.add_nodes_from(posDict)
    G.add_edges_from(semList)
    return(G)

def createEdgeColor(nxGraph):
    """ Create a list of edge colors from the nxGraph edge color attribute. In
    order to get colors working correctly I had to pass nx.draw a list.  """
    edgeColor = [d['color'] for (u,v,d) in nxGraph.edges_iter(data=True)]
    return(edgeColor)

def createEdgeLabel(nxGraph):
    """ Create a list of edge colors from the nxGraph edge color attribute. In
    order to get colors working correctly I had to pass nx.draw a list.  """
    edgeLabel = dict()
    for (u,v,d) in nxGraph.edges_iter(data=True):
        edgeLabel[(u,v)] = d['label']
    return(edgeLabel)

def createNodeColor(nxGraph, expNodes):
    """ Non-expressed isoforms were colored grey to set them apart. """
    nodeColor=list()
    for gene in nxGraph.nodes_iter():
        if gene in expNodes:
            nodeColor.append('white')
        else:
            nodeColor.append('grey')
    return(nodeColor)

if __name__ == "__main__":
    args = getOptions()

    # Turn on Logging if option -g was given
    if args.log:
        setLogger(args.log,logging.INFO)

    # Import Node Positions
    posDict = getNodePos(args.pname)

    # Import SEM
    semList,expNodes = getSemResults(args.sname)

    # Build Network
    mygraph = createGraph(posDict, semList)

    # Pull edge color attributes out so I can use them when graphing
    edgeColor = createEdgeColor(mygraph)

    # Pull edge label attributes out so I can use them when graphing
    edgeLabel = createEdgeLabel(mygraph)

    # Color Nodes that were not expressed grey
    nodeColor = createNodeColor(mygraph, expNodes)
    
    # Create Figure and plot
    logging.info("Plotting")
    fig = plt.figure(figsize=(20,20))
    nx.draw_networkx_nodes(mygraph,pos=posDict,node_size=2500,node_color=nodeColor)
    nx.draw_networkx_edges(mygraph,pos=posDict,edge_color=edgeColor)
    nx.draw_networkx_labels(mygraph,pos=posDict)
    nx.draw_networkx_edge_labels(mygraph,pos=posDict,edge_labels=edgeLabel,label_pos=.5)
    plt.axis('off')
    fig.savefig(args.out, format='png',bbox_inches='tight')

    logging.info("Script Complete")
