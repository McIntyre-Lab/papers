#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:53:16 2026

@author: mcintyre
"""

import pandas as pd
import networkx as nx

#from create_nodes_and_edges_4_genesetid.sas


nodes_df = pd.read_csv("/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/geneset_nodes.csv", dtype={"node_id":str})
edges_df = pd.read_csv("/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/geneset_edges.csv", dtype={"component_id":str,"geneID":str,"source":str})

print(nodes_df.columns)
print(edges_df.columns)

#graph 
G = nx.MultiDiGraph()


# -- add nodes with their attributes -------------------------
for _, row in nodes_df.iterrows():
    G.add_node(row['nodeID'])          

# -- add edges ------------------------------------------------
for _, row in edges_df.iterrows():
    # edge from componentID to geneID
    G.add_edge(row['component_id'],
               row['geneID'],
              
              )  
    
weak_components = list(nx.weakly_connected_components(G))

component_id_map = {}                      # node -> component id
for idx, comp_set in enumerate(weak_components):
    for node in comp_set:
        component_id_map[node] = idx

# Attach component ID as a node attribute
nx.set_node_attributes(G, component_id_map, name='genesetid')
comp_summary = (
    pd.DataFrame({
        'genesetid': [c for c in range(len(weak_components))],
        'num_nodes':    [len(c) for c in weak_components],
                })
                )
print(comp_summary)

geneset_df = pd.DataFrame({
    'nodeID'      : list(G.nodes),
    'genesetID' : [component_id_map[n] for n in G.nodes]
})

#output the link between the nodes and the geneset
geneset_df.to_csv('/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/zenodo/fiveSpecies_network_files/nodes_with_geneset2.csv', index=False)
print("Saved a per‑node CSV: component_map_by_node2.csv")

