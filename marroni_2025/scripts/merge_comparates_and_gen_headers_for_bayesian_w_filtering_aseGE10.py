#!/usr/bin/env python3

import argparse
import os
import sys
import numpy as np
from functools import reduce
from collections import OrderedDict
import pandas as pd

## merge filtered/summarized files with qsim values by user-specified comparison

def getOptions():
    parser = argparse.ArgumentParser(description='Merges together filtered/summarized comparate count tables by user-specified comparison')
    parser.add_argument("-output", "--output", dest="output", action="store", required=True, help="Output directory for complete merged comparate files ready for Bayesian")
    parser.add_argument("-comp", "--comp", dest="comp", action='store', required=True, help="Input filtered/summarized count tables per one comparate")
    parser.add_argument("-design", "--design", dest="design", action='store', required=True, help="Design file")
    args=parser.parse_args()
    return(args)

def main():
    args = getOptions()

    ### Read in design file as dataframe
    df_design = pd.read_csv(args.design)

    ### Create subset of design file of comparate specification columns (will quantify # comparates by number of columns left)
    ### Store compID to name output file

    c1_list = df_design['Comparate_1'].tolist()
    c2_list = df_design['Comparate_2'].tolist()

    dict = {}
    col_list = list(df_design.columns.values)
    row_list = []
    comparison_list = df_design['compID'].tolist()
    del df_design['compID']

    ### Create dictionaries per design file row to store the row's comparate files 
    for index, sample in df_design.iterrows():
        dict[index] = list(sample)

    ## If there are comparison columns (column # > 1)
    for key in dict:
        row_list = dict[key]
        file_list = []
        comp_dict = {}
        comparison = comparison_list[key]
        c1= c1_list[key]
        c2= c2_list[key]
        for i, comp in enumerate(row_list):
            comp_dict[i+1] = comp

            ### Assign filename so it can be called
            row_list[i] = args.comp + '/bayesian_input_' + comp + '.csv'

            file = pd.read_csv(row_list[i], index_col=None, header =0)
            file_list.append(file)


        df_merged = reduce(lambda x, y: pd.merge(x, y, on = ['FEATURE_ID']), file_list)

        ### drop columns you don't want before merge
        df_merged = df_merged[df_merged.columns.drop(list(df_merged.filter(regex='comp')))]

        df_merged.set_index('FEATURE_ID', inplace=True)

        ## AMM fixing below line get_values is deprecated
        ## merged_headers = list(df_merged.columns.get_values())
        merged_headers = list(df_merged.columns.to_numpy())
        print("printing merged headers")
        print(merged_headers)

        ### For stan model, requires headers to have general comparate input names
        ### This reassigns comparate names to be c1, c2, c3... depending on design file specifications
        for x in comp_dict:
            for i in range(len(merged_headers)):
                if c1 in merged_headers[i]:
                   merged_headers[i] = merged_headers[i].replace(c1, 'c1') 
                if c2 in merged_headers[i]:
                   merged_headers[i] = merged_headers[i].replace(c2, 'c2') 

        df_merged.columns=merged_headers


        df_filtered = df_merged

        ## dump genes where c1_flag_analyze and c2_flag_analyze are both 0
        df2_filtered = df_filtered[(df_filtered['c1_flag_analyze']== 1) & (df_filtered['c2_flag_analyze'] == 1)]

        ## dump genes where priors are 0
        df3_filtered = df2_filtered[(df2_filtered['prior_c1_g1']> 0) &
                  (df2_filtered['prior_c1_g2']> 0) &
                  (df2_filtered['prior_c2_g1']> 0) &
                  (df2_filtered['prior_c2_g2']> 0) ].copy()

        ## dump genes where sum of the ase counts for each rep is less than 10
        for comp in ["c1", "c2"]:
            cols = [c for c in df3_filtered.columns if "counts_" + comp + "_g" in c]
            reps = list(set([r.split("_")[-1] for r in cols]))
            for rep in reps:
                sumcols = [c for c in cols if comp in c and rep in c]
                df3_filtered['flag_ge_10_' + comp + '_' + rep ] = np.where(df3_filtered[sumcols].sum(axis=1) >=10, 1, 0)

        gecols = [c for c in df3_filtered.columns if 'flag_ge_10' in c ]
        df3_filtered['flag_ge_10_all'] = np.where(df3_filtered[ gecols ].sum(axis= 1) == len(gecols), 1, 0)

        df4_filtered = df3_filtered[df3_filtered[ gecols ].sum(axis= 1) == len(gecols)].drop(columns=gecols)
       
        outfile = args.output + '/bayesian_input_' + comparison + '.csv'
        df4_filtered.to_csv(outfile)

if __name__=='__main__':
    main()
