import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as Rectangle
import numpy as np
import networkx as nx
import pydot


def importPos(pname):
    """ Import node positions show network structure will match biological pathway """
    mypos=dict()
    with open(pname, 'r') as CSV:
        for row in CSV:
            cols=row.rstrip().split(',')
            mypos[cols[0]] = np.array([float(cols[1]),float(cols[2])])
    return(mypos)


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


def main(pname,dname1,dname2,oname):
    # Pull in node positions from my positions CSV
    mypos = importPos(pname)

    # Create a basic figure and subplot
    fig = plt.figure(figsize=(20,20))
    ax = fig.add_subplot(111)

    # Pull in network structure from dot file
    G1 = nx.read_dot(dname1)

    # Read in factor flags from SAS output
    mycol = importFactor(dname2)

    # Figure out how many factors there are
    numFactors = max(mycol.values())

    # Color pallet based on the jet colors
    jet = plt.cm.jet
    bounds = np.linspace(0,jet.N,numFactors+1)

    # Color each factor differently on the network
    # Add factor color attribute to each node
    colArray = list()
    for node in G1.nodes():
        try:
            colArray.append(bounds[mycol[node]])
        except:
            colArray.append(0)

    # Draw network colors
    nx.draw_networkx(G1,mypos,node_size=3000,font_size=12,node_color=colArray,sytle='solid',cmap='jet',font_color='w',font_weight='bold')

    # Create Legend and Plot
    p, labels = buildLegend(numFactors, jet, bounds)
    plt.legend(p,labels)
    plt.axis('off')
    plt.savefig(oname)


if __name__ == "__main__":
    pname = '/home/jfear/mclab/cegs_sem_sd_paper/documentation/networks/sex_det_node_positions_w_isoforms.csv'
    dname1 = '/home/jfear/mclab/cegs_sem_sd_paper/documentation/networks/sex_det_w_isoforms.dot'
    dname2 = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.csv'
    oname = '/home/jfear/mclab/cegs_sem_sd_paper/analysis_output/factor/dsrp_sex_det_factor_analysis_combined_iso.svg'
    main(pname,dname1,dname2,oname)
