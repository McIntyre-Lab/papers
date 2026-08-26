#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 10:11:54 2025

@author: mcintyre
"""
proj ="/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary/fiveSpecies_annotations/link_files"  # folder that holds the network

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


#nodes_file   => 'node_list.csv'    jxnhash, geneid, source
#edges_file   => 'edges.csv'    source_jxnhash, target_jxnhash, source, target

nodes_df = pd.read_csv("/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/node_list.csv", dtype=str)
edges_df = pd.read_csv("/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/edges.csv", dtype=str)

#build the graph

G = nx.MultiDiGraph()


# add nodes with their attributes 
for _, row in nodes_df.iterrows():
    G.add_node(row['jxnhash'],
               geneid=row['geneID'],
               species=row['source'])          

# add edges 
for _, row in edges_df.iterrows():
    # edge from source_jxnhash to target_jxnhash
    G.add_edge(row['source_jxnHash'],
               row['target_jxnHash'],
               species=row['source'])  # species for the edge (if needed)

#   Identify connected components (ignoring direction)

weak_components = list(nx.weakly_connected_components(G))

component_id_map = {}                      # node -> component id
for idx, comp_set in enumerate(weak_components):
    for node in comp_set:
        component_id_map[node] = idx

# Attach the unique component ID as a node attribute
nx.set_node_attributes(G, component_id_map, name='component_id')
comp_summary = (
    pd.DataFrame({
        'component_id': [c for c in range(len(weak_components))],
        'num_nodes':    [len(c) for c in weak_components],
                })
                )
print(comp_summary)

#  output the link between nodes and components 
nodes_df['component_id'] = nodes_df['jxnhash'].map(component_id_map)
nodes_df.to_csv('nodes_with_components.csv', index=False)

print(f"Nodes: {nodes_df.shape[0]}, Edges: {edges_df.shape[0]}")

node_df = pd.DataFrame({
    'jxnhash'      : list(G.nodes),
    'component_id' : [component_id_map[n] for n in G.nodes]
})

node_df['geneid'] = [G.nodes[n].get('geneid') for n in G.nodes]
node_df['source'] = [G.nodes[n].get('source') for n in G.nodes]

nodes_df.to_csv('/nfshome/mcintyre/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/submission/supplementary_files/fiveSpecies_annotations/link_files/component_map_by_node.csv', index=False)
print("Saved a per‑node CSV: component_map_by_node.csv")

#this is too big to plot nodes=379,602, Edges: 1,419,705

