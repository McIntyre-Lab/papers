#!/usr/bin/env python

import argparse
import pandas as pd
import itertools
from multiprocessing import Pool


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Calculate pairwise distance measures for all pairs of transcripts within the same gene")

    # Input data
    parser.add_argument("-i", "--input-merged", dest="inMerged", required=True, help="CSV file of merged transcript-level junction, fragment, and exon region variables")
    parser.add_argument("-n", "--cpu", dest="inCPU", type=int, default=8, help="Number of CPUs to use for multi-threading of pairwise comparisons (default: 8)")
    # Output data
    parser.add_argument("-o", "--output", dest="outFile", required=True, help="Output CSV file of pairwise transcript distance measures within each gene")

    args = parser.parse_args()
    return args

def calculate_distance(mergeAllDF,cols,genes,n,pairwiseOut=None):
    if n > 1:
        groupResult = pd.DataFrame(columns=cols)
    for gene in genes:
        # Get list of all transcripts in the gene
        xcrptList = list(mergeAllDF.loc[mergeAllDF['gene_id']==gene,'transcript_id'])
        # Skip single transcript genes
        if len(xcrptList) == 1:
            continue
        # Generate list of pairs
        pairList = list(itertools.combinations(xcrptList,2))
        # Calculate distances for each pair in gene (genes are parallelized using multiprocessing)
        for T1,T2 in pairList:
            t1DF = mergeAllDF.loc[mergeAllDF['transcript_id']==T1].squeeze()
            t2DF = mergeAllDF.loc[mergeAllDF['transcript_id']==T2].squeeze()
            # Make empty Series to place distance values
            singlePair = pd.Series(index=cols)
            singlePair[['gene_id','transcript_1','transcript_2']] = gene,T1,T2
            
            ### Get juction distance values
            # Both transcripts are monoexon
            if t1DF['flag_monoexon_transcript'] == 1 and t2DF['flag_monoexon_transcript'] == 1:
                singlePair[['num_junction_T1_only','num_junction_T2_only','num_junction_shared']] = 0
                singlePair['prop_junction_diff'] = 0
                singlePair['prop_junction_similar'] = 1
                singlePair[['junction_T1_only','junction_T2_only','junction_shared']] = ""
            # Only T1 is monoexon
            elif t1DF['flag_monoexon_transcript'] == 1 and t2DF['flag_monoexon_transcript'] == 0:
                singlePair[['num_junction_T1_only','num_junction_shared']] = 0
                singlePair['prop_junction_similar'] = 0
                singlePair['prop_junction_diff'] = 1
                singlePair[['junction_T1_only','junction_shared']] = ""
                singlePair['num_junction_T2_only'] = len(t2DF['junctionID_order'].split("|"))
                singlePair['junction_T2_only'] = t2DF['junctionID_order']
            # Only T2 is monoexon
            elif t1DF['flag_monoexon_transcript'] == 0 and t2DF['flag_monoexon_transcript'] == 1:
                singlePair[['num_junction_T2_only','num_junction_shared']] = 0
                singlePair['prop_junction_similar'] = 0
                singlePair['prop_junction_diff'] = 1
                singlePair[['junction_T2_only','junction_shared']] = ""
                singlePair['num_junction_T1_only'] = len(t1DF['junctionID_order'].split("|"))
                singlePair['junction_T1_only'] = t1DF['junctionID_order']
            # Both transcripts multi-exon
            else:
                t1Junc = t1DF['junctionID_order'].split("|")
                t2Junc = t2DF['junctionID_order'].split("|")
                juncSharedSet = set(t1Junc).intersection(set(t2Junc))
                juncT1Set = set(t1Junc).difference(set(t2Junc))
                juncT2Set = set(t2Junc).difference(set(t1Junc))
                singlePair['num_junction_T1_only'] = len(juncT1Set)
                singlePair['num_junction_T2_only'] = len(juncT2Set)
                singlePair['num_junction_shared'] = len(juncSharedSet)
                singlePair['prop_junction_diff'] = (singlePair['num_junction_T1_only'] + singlePair['num_junction_T2_only'])/(singlePair['num_junction_T1_only'] + singlePair['num_junction_T2_only'] + singlePair['num_junction_shared'])
                singlePair['prop_junction_similar'] = 1 - singlePair['prop_junction_diff']
                singlePair['junction_T1_only'] = "|".join(sorted(juncT1Set, key=t1Junc.index))
                singlePair['junction_T2_only'] = "|".join(sorted(juncT2Set, key=t2Junc.index))
                singlePair['junction_shared'] = "|".join(sorted(juncSharedSet, key=t2Junc.index))
        
            ### Get exon region distance values
            t1ER = t1DF['exonRegionID_order'].split("|")
            t2ER = t2DF['exonRegionID_order'].split("|")
            ERSharedSet = set(t1ER).intersection(set(t2ER))
            ERT1Set = set(t1ER).difference(set(t2ER))
            ERT2Set = set(t2ER).difference(set(t1ER))
            singlePair['num_ER_T1_only'] = len(ERT1Set)
            singlePair['num_ER_T2_only'] = len(ERT2Set)
            singlePair['num_ER_shared'] = len(ERSharedSet)
            singlePair['prop_ER_diff'] = (singlePair['num_ER_T1_only'] + singlePair['num_ER_T2_only'])/(singlePair['num_ER_T1_only'] + singlePair['num_ER_T2_only'] + singlePair['num_ER_shared'])
            singlePair['prop_ER_similar'] = 1 - singlePair['prop_ER_diff']
            singlePair['ER_T1_only'] = "|".join(sorted(ERT1Set, key=t1ER.index))
            singlePair['ER_T2_only'] = "|".join(sorted(ERT2Set, key=t2ER.index))
            singlePair['ER_shared'] = "|".join(sorted(ERSharedSet, key=t2ER.index))
        
            ### Get fragment distance values
            t1Frag = t1DF['fragmentID_order'].split("|")
            t2Frag = t2DF['fragmentID_order'].split("|")
            fragSharedSet = set(t1Frag).intersection(set(t2Frag))
            fragSharedSingSet = set(x for x in fragSharedSet if x.startswith('S'))
            fragT1Set = set(t1Frag).difference(set(t2Frag))
            fragT1SingSet = set(x for x in fragT1Set if x.startswith('S'))
            fragT2Set = set(t2Frag).difference(set(t1Frag))
            fragT2SingSet = set(x for x in fragT2Set if x.startswith('S'))
            singlePair['num_fragment_T1_only'] = len(fragT1Set)
            singlePair['num_fragment_T2_only'] = len(fragT2Set)
            singlePair['num_fragment_shared'] = len(fragSharedSet)
            singlePair['prop_fragment_diff'] = (singlePair['num_fragment_T1_only'] + singlePair['num_fragment_T2_only'])/(singlePair['num_fragment_T1_only'] + singlePair['num_fragment_T2_only'] + singlePair['num_fragment_shared'])
            singlePair['prop_fragment_similar'] = 1 - singlePair['prop_fragment_diff']
            singlePair['fragment_T1_only'] = "|".join(sorted(fragT1Set, key=t1Frag.index))
            singlePair['fragment_T2_only'] = "|".join(sorted(fragT2Set, key=t2Frag.index))
            singlePair['fragment_shared'] = "|".join(sorted(fragSharedSet, key=t2Frag.index))
            singlePair['num_fragment_singletons_T1_only'] = len(fragT1SingSet)
            singlePair['num_fragment_singletons_T2_only'] = len(fragT2SingSet)
            singlePair['num_fragment_singletons_shared'] = len(fragSharedSingSet)
            
            ### Get nt distance values
            t1FragNT = list(map(int,t1DF['fragmentLength_order'].split("|")))
            t2FragNT = list(map(int,t2DF['fragmentLength_order'].split("|")))
            singlePair['num_nt_shared'] = sum([int(t1FragNT[t1Frag.index(x)]) for x in t1Frag if x in fragSharedSet])
            singlePair['num_nt_T1_only'] = sum([int(t1FragNT[t1Frag.index(x)]) for x in t1Frag if x in fragT1Set])
            singlePair['num_nt_T2_only'] = sum([int(t2FragNT[t2Frag.index(x)]) for x in t2Frag if x in fragT2Set])
            singlePair['total_nt'] = singlePair['num_nt_shared'] + singlePair['num_nt_T1_only'] + singlePair['num_nt_T2_only']
            singlePair['prop_nt_diff'] = (singlePair['num_nt_T1_only'] + singlePair['num_nt_T2_only'])/(singlePair['total_nt'])
            singlePair['prop_nt_similar'] = 1 - singlePair['prop_nt_diff']
            singlePair['num_nt_T1_only_in_shared_ER'] = sum([int(t1FragNT[t1Frag.index(x)]) for x in t1Frag if (x in fragT1Set) and (x.split(":")[0] in ERSharedSet)])
            singlePair['num_nt_T2_only_in_shared_ER'] = sum([int(t2FragNT[t2Frag.index(x)]) for x in t2Frag if (x in fragT2Set) and (x.split(":")[0] in ERSharedSet)])
            singlePair['num_nt_shared_in_shared_ER'] = sum([int(t2FragNT[t2Frag.index(x)]) for x in t2Frag if (x in fragSharedSet) and (x.split(":")[0] in ERSharedSet)])
            singlePair['total_nt_in_shared_ER'] = singlePair['num_nt_shared_in_shared_ER'] + singlePair['num_nt_T1_only_in_shared_ER'] + singlePair['num_nt_T2_only_in_shared_ER']
            if singlePair['total_nt_in_shared_ER'] != 0:
                singlePair['prop_nt_diff_in_shared_ER'] = (singlePair['num_nt_T1_only_in_shared_ER'] + singlePair['num_nt_T2_only_in_shared_ER'])/(singlePair['total_nt_in_shared_ER'])
                singlePair['prop_nt_similar_in_shared_ER'] = 1 - singlePair['prop_nt_diff_in_shared_ER']
            else:
                singlePair['prop_nt_diff_in_shared_ER'] = 0
                singlePair['prop_nt_similar_in_shared_ER'] = 0
            singlePair['num_nt_T1_only_in_unique_ER'] = sum([int(t1FragNT[t1Frag.index(x)]) for x in t1Frag if (x in fragT1Set) and (x.split(":")[0] in ERT1Set)])
            singlePair['num_nt_T2_only_in_unique_ER'] = sum([int(t2FragNT[t2Frag.index(x)]) for x in t2Frag if (x in fragT2Set) and (x.split(":")[0] in ERT2Set)])

            ### Append distance measures of pair gene result
            if n == 1:
                pairwiseOut.write(singlePair.to_frame().T.to_csv(header=False, index=False))
            else:
                groupResult = groupResult.append(singlePair.to_frame().T, ignore_index=True)
    if n > 1:
        return groupResult


def chunks(lst, n):
    # Split list into n chunks and return list of lists
    splitSize = (len(lst)//n) + ((len(lst)%n)>0)
    list_of_lists = []
    for element in range(0, len(lst), splitSize):
        list_of_lists.append(lst[element:element + splitSize])
    return list_of_lists

result_list = []
def callback_results(result):
    # Callback function to append result to list of results
    result_list.append(result)
    
def main():
    # Get input merged files
    mergeAllDF = pd.read_csv(args.inMerged, low_memory=False)
    
    cols=['gene_id','transcript_1','transcript_2','num_junction_T1_only','num_junction_T2_only',
             'num_junction_shared','prop_junction_diff','prop_junction_similar',
             'junction_T1_only','junction_T2_only','junction_shared','num_ER_T1_only',
             'num_ER_T2_only','num_ER_shared','prop_ER_diff','prop_ER_similar',
             'ER_T1_only','ER_T2_only','ER_shared','num_fragment_T1_only',
             'num_fragment_T2_only','num_fragment_shared','prop_fragment_diff',
             'prop_fragment_similar','fragment_T1_only','fragment_T2_only','fragment_shared',
             'num_fragment_singletons_T1_only','num_fragment_singletons_T2_only',
             'num_fragment_singletons_shared','num_nt_shared','num_nt_T1_only',
             'num_nt_T2_only','total_nt','prop_nt_diff','prop_nt_similar',
             'num_nt_T1_only_in_shared_ER','num_nt_T2_only_in_shared_ER',
             'num_nt_shared_in_shared_ER','total_nt_in_shared_ER','prop_nt_diff_in_shared_ER',
             'prop_nt_similar_in_shared_ER','num_nt_T1_only_in_unique_ER',
             'num_nt_T2_only_in_unique_ER']

    # If 1 cpu available
    if args.inCPU == 1:
        pairwiseOut = open(args.outFile, 'w')
        pairwiseOut.write(",".join(cols)+"\n")
        # Loop through all genes and calculate distances
        calculate_distance(mergeAllDF,cols,mergeAllDF['gene_id'].unique(),args.inCPU,pairwiseOut)

    # If cpu > 1, parallelize
    elif args.inCPU > 1:
        # Get lists for each process based on cpu value
        geneLists = chunks(mergeAllDF['gene_id'].unique(),args.inCPU)
        # Generate multiprocess Pool with specified number of cpus
        #     to loop through genes and calculate distances
        pool = Pool(args.inCPU)
        for genes in geneLists:
            geneInfo = mergeAllDF[mergeAllDF['gene_id'].isin(genes)]
            pool.apply_async(calculate_distance, args=(geneInfo,cols,genes,args.inCPU), callback=callback_results)
        pool.close()
        pool.join()        
        # Output results
        pd.concat(result_list).to_csv(args.outFile, index=False)
    else:
        print("ERROR: Invalid cpu parameter")

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()

