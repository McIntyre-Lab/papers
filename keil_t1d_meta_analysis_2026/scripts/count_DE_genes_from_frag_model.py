#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


@author: nkeil
"""

import pandas as pd
#import seaborn as sns
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse


def getOptions():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Analyze DE Gene Counts')
    parser.add_argument('--in_file', type=str, required=True, help='Path to the input file from SAS model')
    parser.add_argument('--effect', type=str, required=True, help='Effect for which the model is being checked for DE')
    parser.add_argument('--t1d', type=str, required=True, help='Path to the T1D gene list')
    parser.add_argument('--immune', type=str, required=True, help='Path to the immune gene list')
    parser.add_argument('--pvalue', type=float, default=0.001, help='pvalue for determining if a fragment is significant')
    parser.add_argument('--prop', type=float, default=0.2, help='Proportion cutoff for flagging DE genes')
    parser.add_argument('--out_graph', type=str, required=True, help='Filename for saving the output graphs')
    parser.add_argument('--out_csv', type=str, required=True, help='Filename for saving the output CSV of DE gene counts')

    args = parser.parse_args()
    return args


def main():
    # Read the input files
    t1DF = pd.read_csv(args.in_file, sep=',', header=0, low_memory=False)
    t1DF['geneID'] = t1DF['featureID'].str.split(':').str[0]
    
    #Dealing with <0.0000000001
    t1DF['ProbF'] = pd.to_numeric(t1DF['ProbF'].str.replace('<', ''), errors='coerce')
    
    ##Flag significant fragments based on pvalue input
    sig=args.pvalue
    
    t1DF['flag_sig'] = (t1DF['ProbF'] < sig).astype(int)
    
    #Count significant fragments
    gene_effect_DF = t1DF.groupby(['geneID', 'Effect']).agg(
        cnt_frags_sig=('flag_sig', 'sum'),
        total_frags=('geneID', 'size')
    )
    
    #Gte proportion of frgaments significant
    gene_effect_DF['prop_frags_sig'] = gene_effect_DF['cnt_frags_sig'] / gene_effect_DF['total_frags']
    gene_effect_DF.reset_index(inplace=True)
    
    #Subset DF based on effect of interest
    effect_DF = gene_effect_DF[gene_effect_DF['Effect'] == args.effect]
    
    # Flag for proportion_of_frags_sig greater than prop
    effect_DF['flag_gene_DE'] = (effect_DF['prop_frags_sig'] > args.prop).astype(int)
    
    # Import T1D and immune gene lists
    t1d_geneDF = pd.read_csv(args.t1d, sep=',', header=0, low_memory=False)
    immune_geneDF = pd.read_csv(args.immune, sep=',', header=0, low_memory=False)
    
    effect_DF['flag_T1D'] = effect_DF['geneID'].isin(t1d_geneDF['ensembl_geneID']).astype(int)
    effect_DF['flag_immune'] = effect_DF['geneID'].isin(immune_geneDF['Ensembl_ID']).astype(int)
    
    # Calculate DE gene counts
    cnt_DE_all = effect_DF['flag_gene_DE'].sum()
    cnt_DE_immune = (effect_DF['flag_gene_DE'] * effect_DF['flag_immune']).sum()
    cnt_DE_t1d = (effect_DF['flag_gene_DE'] * effect_DF['flag_T1D']).sum()
    
    #Get counts of analyzable genes
    cnt_all = len(effect_DF)
    cnt_immune = effect_DF['flag_immune'].sum()
    cnt_T1D = effect_DF['flag_T1D'].sum()
    
    
    #Calculate proportions
    prop_DE_all = cnt_DE_all / cnt_all
    prop_DE_immune = cnt_DE_immune / cnt_immune if cnt_immune > 0 else 0
    prop_DE_t1d = cnt_DE_t1d / cnt_T1D if cnt_T1D > 0 else 0
    
    #Proportions data for plotting
    proportions = [prop_DE_all, prop_DE_immune, prop_DE_t1d]
    categories = ['All Genes', 'Immune Genes', 'T1D Genes']
    
    #Plot proportions
    
    plt.figure(figsize=(7, 6))
    plt.style.use('ggplot')
    
    plt.bar(categories, proportions, color='skyblue')
    plt.title('Proportions of DE Genes - ' + args.effect + " effect")
    plt.ylabel('Proportion of DE Genes')
    plt.ylim(0, 1)
    
    #Save the plot
    plt.tight_layout()
    plt.savefig(args.out_graph)
    
    #Create a DataFrame for DE gene counts
    de_gene_counts = pd.DataFrame({
    'Category': ['All', 'Immune', 'T1D'],
    'DE_gene_count': [cnt_DE_all, cnt_DE_immune, cnt_DE_t1d],
    'Prop_genes_DE': [prop_DE_all, prop_DE_immune, prop_DE_t1d],
    })
    
    #Save the DE gene counts to a CSV file
    de_gene_counts.to_csv(args.out_csv, index=False)
    
    #Output analyzable gene counts to log
    print("Analyzable genes")
    print("All genes: " + str(cnt_all))
    print("Immune genes: " + str(cnt_immune))
    print("T1D genes: " + str(cnt_T1D))
    
if __name__ == '__main__':
    #Parse command line arguments
    global args
    args = getOptions()
    main()