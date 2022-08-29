#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 14:46:34 2021

@author: zach
"""

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from upsetplot import UpSet


def direction_determine(x):
    for i in x[1:]:
        if x[0] * i <=0:
            return 0
    return 1
        

def all_sig(x):
    for i in x:
        if i != 1:
            return 0
    return 1

#for group in ['TCA', 'UGT', 'NI', 'NInoN2']:
#for group in ['NInoN2']:
for group in ['ALL']:

    if group == 'TCA':
        muts = ["KJ550", "RB2347", "VC1265", "AUM2073", "VC2524"]
    elif group == 'UGT':
        muts = ["RB2055", "RB2550", "RB2011", 'UGT49', 'UGT60']
    elif group == 'NI':
        muts = ["CB4856", "CX11314", "DL238", "N2"]
    elif group == 'NInoN2':
        muts = ["CB4856", "CX11314", "DL238"]
    elif group == 'COMBO':
        muts = ["VC1265", "DL238", "UGT60"]
    elif group == 'ALL':
        muts = ["AUM2073","CB4856","CX11314","DL238","KJ550","N2","RB2011","RB2055","RB2347","RB2550","UGT49","UGT60","VC1265","VC2524"]
    else:
        print("ERROR: Invalid group. Must be UGT, TCA, NI or NInoN2.")

    rank = 'rank'
    for type in ['rp_pos', 'rp_neg', 'hilic_pos']:

        ## by pathway
        combined = pd.read_csv("~/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_analysis_" + type + "/MA_FE_rank_byPath_" + type + "_" + group + "_SMD_summary.tsv", sep='\t')
        combined.columns = ["uniqueID"] + combined.columns[1:].tolist()
        combined = combined[["uniqueID", 'effect', 'p_value']]
        combined.columns = ['uniqueID', 'effect_meta_pathway', 'p_meta_pathway']

        df_mut = combined
        mutpath = "~/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_analysis_" + type + "/"
        ## by mutant
        for mut in muts:
            f = mutpath + "MA_FE_rank_byMut_" + type + "_" + mut + "_SMD_summary.tsv"
            dat = pd.read_csv(f)
            dat.columns =  ["uniqueID"] + dat.columns[1:].tolist()
            dat = dat[['uniqueID', 'effect', 'p_value']]
            column_names = ['uniqueID'] + ['effect_' + mut] + ['p_' + mut]
            dat.columns = column_names
            df_mut = df_mut.merge(dat, on = 'uniqueID', how = 'left')

        df_mut.index = df_mut["uniqueID"]
        df_mut = df_mut.drop("uniqueID", axis = 1)

        es_cols = pd.Series(df_mut.columns).str.contains("effect|uniqueID").tolist()
        df_es = df_mut.loc[:,es_cols]

        p_cols = pd.Series(df_mut.columns).str.contains("p_|uniqueID").tolist()
        df_p = df_mut.loc[:,p_cols]


        df_es["flag_same_direction"] = df_es.apply(direction_determine, axis = 1)

        ##### result 2: counts same direction
        counts_same_direct = sum(df_es["flag_same_direction"])
 
    #####   
        sig_05_table = pd.DataFrame(np.where(df_p < 0.05, 1, 0))
        sig_05_table.columns = df_p.columns
        sig_05_table.index = df_p.index


    ##### result 1: count table for sig 0.05
        counts_sig_05 = pd.DataFrame(sig_05_table.apply(sum, axis = 0))
        counts_sig_05.columns = ['c_sig_0.05']
 
    #####
        both_sig_same_05 = pd.Series(sig_05_table.loc[(df_es["flag_same_direction"] == 1) & (sig_05_table['p_meta_pathway'] == 1)].index)
        both_sig_same_05 = sig_05_table.apply(lambda x: len(df_p.loc[(df_es["flag_same_direction"] == 1) & (x ==1), x.name]))
        both_sig_same_05.name = 'both_sig_0.05_and_same_direct'

    #####

    ####

        sig_05_table_combine = sig_05_table[sig_05_table['p_meta_pathway'] == 1]
  

        df_es_sig_05_combine = df_es[sig_05_table["p_meta_pathway"] == 1]


        df_es_sig_05_combine_sdirect = df_es_sig_05_combine[df_es_sig_05_combine["flag_same_direction"]==1]
        print(len(df_es_sig_05_combine_sdirect))
        df_es_sig_05_combine_ndirect = df_es_sig_05_combine[df_es_sig_05_combine["flag_same_direction"]==0][100:200]
#######################
############ extract a effect size matrix for heatmap in R
        es_matrix = df_es[df_es['flag_same_direction'] == 1]
        es_matrix['sig_counts'] = sig_05_table.loc[df_es['flag_same_direction'] == 1, ~sig_05_table.columns.str.contains("meta_pathway")].apply(sum, axis = 1)
        es_matrix['direction'] = np.where(es_matrix['effect_meta_pathway'] > 0, 1, 0)
        es_matrix = es_matrix.sort_values(by = ['sig_counts', 'direction'], ascending = False)
        sort_info = es_matrix[['sig_counts', 'direction']]
        #es_matrix = es_matrix.drop(['flag_same_direction', 'sig_counts', 'direction'], axis = 1)
        es_matrix = es_matrix.merge(df_p["p_meta_pathway"], left_index = True, right_index = True, how = 'left')
        es_matrix.to_csv("~/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/effect_size_matrix_hmp_" + type + "_" + group + ".csv")
        print(len(es_matrix))



####### scattor all same direction
        fig = plt.figure(figsize=(15,4), facecolor='white',edgecolor='black')
        ax = fig.add_subplot(111)
#        colors = ['red', 'blue', 'green', 'orange', 'black', 'purple']
        colors = sns.color_palette("tab20", 20).as_hex()
        print(colors)

        labels = df_es_sig_05_combine.columns[:-1]
        labels = pd.Series(labels).str.split("_").apply(lambda x: x[1])
        sizes = [4] + [2] * 19
        s = 0
        patches = []
        for i in df_es_sig_05_combine.columns[:-1]:
            patch = ax.plot([i * 2 for i in range(0, len(df_es_sig_05_combine_sdirect))], 
                df_es_sig_05_combine_sdirect.loc[:, i].tolist(),
                'o', markersize = sizes[s], color = colors[s])
            s +=1
            patches += patch
        ax.hlines(0, ax.get_xlim()[0], ax.get_xlim()[1], color = 'black', linewidth = 0.5)
        ax.legend([matplotlib.lines.Line2D([0],[0], color = x, marker = 'o', ls = '', markersize=3) for x in colors], 
                  labels, loc = 1, prop = dict(size = 6))
        fig.savefig('/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/effect_all_model_same_direct_' + type + "_" + group + '.pdf')
        plt.close(fig)

    ##### upset plot
        def plot_upset(df,title, fig, boxCols=None):
            """
            Plot an UpSet plot given a properly formatted dataframe (multi-indexed with boolean values)
            """
            upset = UpSet(df,subset_size='count',show_counts=True,sort_by='degree',sort_categories_by=None)
            if boxCols is not None:
                if type(boxCols) != list:
                    boxCols = list(boxCols)
                for col in boxCols:
                    upset.add_catplot(value=col, kind='box', elements = 4, showfliers=False, color=sns.color_palette('colorblind',15).as_hex()[boxCols.index(col)])
            upset.plot(fig = fig)
            plt.subplots_adjust(right=1.00001)
            plt.suptitle(title)
    
        def split_pos_neg(sig_table, es_table):
            cols = sig_table.columns
            res = pd.DataFrame(index = sig_table.index)
            for col in cols:
                sig_list = sig_table[col]
                es_list = es_table[col]
                res[col + "_pos"] = 0  ## _ge_PD1074
                res[col + "_neg"] = 0  ## _le_PD1074
           
                res[col+"_pos"] = np.where(sig_list & es_list ,1,0) + res[col+'_pos']
                res[col+"_neg"] = np.where((sig_list == 1) & (es_list == 0), 1, 0) + res[col+'_neg']
            return res
        
    
    

### 
        sig_upset = sig_05_table.iloc[:,1:]
        sig_upset = sig_upset[~(df_es['flag_same_direction'] == 1)]
        sig_upset = sig_upset[~(sig_upset.apply(sum, axis = 1) == 0)]
###
        df_es_direct = df_es.loc[sig_upset.index, df_es.columns[1:-1]]
        df_es_direct = pd.DataFrame(np.where(df_es_direct > 0, 1, 0))
        df_es_direct.index = sig_upset.index
        df_es_direct.columns = sig_upset.columns
        inputs = split_pos_neg(sig_upset, df_es_direct)

###
        cols = inputs.columns.to_list()
        inputs = inputs.reset_index()
        inputs[cols] = inputs[cols].astype(bool)
        inputs= inputs.set_index(cols)

###
        fig  = plt.figure(figsize=(12,8))
        plot_upset(inputs,"", fig)
        fig.savefig("/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/upset_sigCounts_not_all_same_direct_" + type + "_" + group + ".pdf", format='pdf')
        fig.savefig("/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/upset_sigCounts_not_all_same_direct_" + type + "_" + group + ".svg", format='svg')
        plt.close(fig)
