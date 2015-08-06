import argparse
import logging
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rectangle
import itertools
import numpy as np
import networkx as nx
import pydot

def getOptions():
    """ Function to pull in arguments """
    parser = argparse.ArgumentParser(description="Function to plot Factors and color them.")
    parser.add_argument("--factor", dest="fname", action='store', required=True, help="Formated output from proc factor [Required]")
    parser.add_argument("-o", dest="oname", action='store', required=True, help="Output PNG name [Required]")
    parser.add_argument("-fs", dest="font", action='store', required=False, default=11, type=int, help="Font size [Default 11pt]")
    parser.add_argument("-cs", dest="canvas", action='store', required=False, default=20, type=int, help="Canvas size [Default 20]")
    args = parser.parse_args()
    #args = parser.parse_args(['--factor', '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_all_iso.csv', '-o', '/home/jfear/tmp/test.png'])
    #args = parser.parse_args(['--factor', '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.csv', '-o', '/home/jfear/tmp/test.png'])
    return(args)

def buildPos(names,args):
    """ Function to construct a grid of equally spaced number for plotting """
    # How many genes per row?
    numPerRow = np.ceil(np.sqrt(len(names)))

    # Space everything out equally in a grid
    coords = np.linspace(0, args.canvas, num=numPerRow)
    X,Y = np.meshgrid(coords, coords[::-1])

    # Create dictionary where key is a name and value are X,Y coords
    myDict = dict()
    for index, name in enumerate(names):
        myDict[name] = [X.ravel()[index], Y.ravel()[index]]

    return(myDict)

def importFactor(dname2):
    """ Import factor flags for each gene/isoform """
    mycol = dict()
    with open(dname2,'r') as FF:
        next(FF)
        for row in FF:
            cols = row.rstrip().split(',')
            mycol[cols[0]] = int(cols[1])
    return(mycol)

def buildLegend(numFactors, jet, bounds):
    """ Use rectangles and factor color pallet to build a legend """
    # I pieced together code on the internet. 
    # Build a bunch of rectangles for display in the legend.
    p = list()
    for i in xrange(numFactors+1):
        rect = Rectangle((0,0),1,1,facecolor=jet(int(bounds[i])))
        p.append(rect)

    # Create labels for the legend.
    labels=['None']
    for i in xrange(numFactors):
        labels.append('Factor' + str(i+1))

    return(p, labels)

if __name__ == "__main__":
    # Pull in command line arguments
    args = getOptions()

    # Read in factor flags from SAS output
    myCol = importFactor(args.fname)
    names = myCol.keys()
    names.sort()

    # Figure out positioning
    myPos = buildPos(names, args)

    # Initiate Network and add nodes
    G = nx.Graph()
    G.add_nodes_from(names)

    # Initialize figure and subplot
    fig = plt.figure(figsize=(args.canvas,args.canvas))
    ax = fig.add_subplot(111)

    # Figure out how many factors there are
    numFactors = max(myCol.values())

    # Color pallet 
    colors = plt.cm.Dark2
    bounds = np.linspace(0,colors.N,numFactors+1)

    # Color each factor differently on the network
    # Add factor color attribute to each node
    colArray = list()
    for node in G.nodes():
        try:
            colArray.append(bounds[myCol[node]])
        except:
            colArray.append(0)

    # Draw network colors
    nx.draw_networkx(G,myPos, node_size=3000,font_size=args.font,node_color=colArray,sytle='solid',cmap='Dark2',font_color='k',font_weight='bold')

    # Create Legend and Plot
    p, labels = buildLegend(numFactors, colors, bounds)
    plt.legend(p,labels)
    plt.axis('off')
    plt.savefig(args.oname)
