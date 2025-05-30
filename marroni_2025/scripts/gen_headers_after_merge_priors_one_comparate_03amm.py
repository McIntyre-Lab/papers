import argparse
import os
import pandas as pd
import numpy as np
from functools import reduce

DEBUG = False

def getOptions():
    parser = argparse.ArgumentParser(
        description='Create headers for model for a single comparate.'
    )
    parser.add_argument(
        "-output",
        "--output",
        dest="output",
        action="store",
        required=True,
        help="Output directory for complete merged comparate files ready for Bayesian"
    )
    parser.add_argument(
        "-input_dir",
        "--input_dir",
        dest="input_dir",
        action="store",
        required=True,
        help="Directory containing ASE count table files"
    )
    parser.add_argument(
        "-design",
        "--design",
        dest="design",
        action='store',
        required=True,
        help="Design file"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        default=False,
        help="Print debugging output"
    )
    args = parser.parse_args()
    return args

def main():
    args = getOptions()
    global DEBUG
    if args.debug:
        DEBUG = True

    if DEBUG:
        print(f"DEBUG: Input directory:\n{args.input_dir}")

    ### Read in design file as dataframe
    # Make sure design file is read as a tsv
    df_design = pd.read_csv(args.design, sep=',')
    if DEBUG:
        print(f"DEBUG: design:\n{df_design}")

    ### Subset design file and create comparate specification columns 
    ### Store compID to name output file
    c1_list = df_design['Comparate_1'].tolist()

    sample_dict = {}
    row_list = []
    comparison_list = df_design['compID'].tolist()
    del df_design['compID']

    ### Create dictionaries per design file row to store the row's comparate files 
    for index, sample in df_design.iterrows():
        sample_dict[index] = list(sample)

    ## Create a dictionary containing the comparisons between each parental genome the comparate
    for key in sample_dict:
        row_list = sample_dict[key]
        file_list = []
        comp_dict = {}
        comparison = comparison_list[key]
        c1 = c1_list[key]
        for i, comp in enumerate(row_list):
            comp_dict[i + 1] = comp

            ### Assign filename so it can be called
            row_list[i] = args.input_dir + '/bayesian_input_' + comp + '.csv'

            file = pd.read_csv(row_list[i], index_col=None, header =0)
            file_list.append(file)

            # Construct the file path
#            file_path = os.path.join(args.input_dir, f'{comp}.csv')

            # Check if file exists
#            if not os.path.exists(file_path):
#                print(f"Error: Input file {file_path} does not exist.")
#                exit(1)

            # Use pd.read_table to read file into dataframe
#            file = pd.read_table(file_path, index_col=None, header=0)
#            file_list.append(file)

        df_merged = reduce(lambda x, y: pd.merge(x, y, on=['FEATURE_ID']), file_list)

        ### Drop columns you don't want before merge
        df_merged = df_merged[df_merged.columns.drop(list(df_merged.filter(regex='comp')))]

        df_merged.set_index('FEATURE_ID', inplace=True)

        # Update to use the newer pandas method
        merged_headers = list(df_merged.columns.to_numpy())

        ### For stan model, requires headers to have general comparate input names
        ### This reassigns comparate names to be c1
        for x in comp_dict:
            for i in range(len(merged_headers)):
                if c1 in merged_headers[i]:
                   merged_headers[i] = merged_headers[i].replace(c1, 'c1')

        df_merged.columns = merged_headers

        df_filtered = df_merged

        ## AMM adding filtering below
        ## dump genes where c1_flag_analyze is  0
        df2_filtered = df_filtered[(df_filtered['c1_flag_analyze']== 1)]

        ## dump genes where priors are 0
        df3_filtered = df2_filtered[(df2_filtered['prior_c1_g1']> 0) &
                  (df2_filtered['prior_c1_g2']> 0) ].copy()

        ## dump genes where sum of the counts for each rep is less than 10
        for comp in ["c1"]:
            cols = [c for c in df3_filtered.columns if "counts_" + comp + "_g" in c]
            reps = list(set([r.split("_")[-1] for r in cols]))
            for rep in reps:
                sumcols = [c for c in cols if comp in c and rep in c]
                df3_filtered['flag_ge_10_' + comp + '_' + rep ] = np.where(df3_filtered[sumcols].sum(axis=1) >=10, 1, 0)

        gecols = [c for c in df3_filtered.columns if 'flag_ge_10' in c ]
        df3_filtered['flag_ge_10_all'] = np.where(df3_filtered[ gecols ].sum(axis= 1) == len(gecols), 1, 0)

        df4_filtered = df3_filtered[df3_filtered[ gecols ].sum(axis= 1) == len(gecols)].drop(columns=gecols)


        # Change the output name from bayesian_input_comp to bayesian_input
#        outfile = os.path.join(args.output, f'bayesian_input_{comparison}')
       	outfile = args.output + '/bayesian_input_' + comparison + '_1cond.csv'
        df4_filtered.to_csv(outfile)

if __name__ == '__main__':
    main()

