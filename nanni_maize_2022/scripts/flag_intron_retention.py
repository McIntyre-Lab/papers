#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import gffutils

def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Flag exonic regions (fusions) with intron retention events")

    # Input data
    parser.add_argument("-d", "--db", dest="inDB", required=True, help="Event Analysis converted GFF3 DB file (*.gff3.db)")
    parser.add_argument("-f", "--fragment-annotation", dest="inFrag", required=True, help="Exon fragment annotation CSV file from Event Analysis (*_exon_fragment_annotations.csv)")

    # Output data
    parser.add_argument("-o", "--output-directory", dest="outDir", required=True, help="Output directory")
    parser.add_argument("-p", "--output-prefix", dest="outPrefix", required=True, help="Prefix for output files")

    args = parser.parse_args()
    return args

def split_column_by_sep(df,col_name=None,sep=None,sort_list=None):
    # Split variable by some character like '|' or ',' and keep all other values the same
    if col_name == None:
        col_name = 'transcript_id'
    if sep == None:
        sep = "|"
    splitList = df[col_name].str.split(sep).apply(pd.Series, 1).stack()
    splitList.index = splitList.index.droplevel(-1)
    tempDF = df.copy()
    del(tempDF[col_name])
    splitDF = tempDF.join(splitList.rename(col_name))
    if sort_list != None:
        splitDF = splitDF.sort_values(by=sort_list)
    del(tempDF, splitList)
    return splitDF

def transform_intron_id(feature):
    # Set intron ID value in attributes as concatenation of exonID values with "_" between
    feature.attributes['ID'] = ["_".join(feature.attributes['Name'])]
    return feature

def main():
    # Get gffutils GFF database file
    db = gffutils.FeatureDB(args.inDB)
    
    # Get all introns and make database in memory
    introns = db.create_introns()
    intronDB = gffutils.create_db(introns,":memory:",disable_infer_genes=True,
                                  disable_infer_transcripts=True,force=True,
                                  id_spec='ID',transform=transform_intron_id,
                                  merge_strategy='merge')
    
    # Output list of exons that overlap entirely with introns
    exonDF = pd.DataFrame(columns=['chr','exonID','exon_attribute','exon_start',
                                   'exon_end','intronID','intron_start',
                                   'intron_end','intron_attribute','flag_intron_retention'])
    for exon in db.features_of_type('exon'):
        hasOverlap = False
        for intron in intronDB.region(exon,completely_within=True):
            hasOverlap = True
            exonDF = exonDF.append({'chr':exon.chrom,
                                    'exonID':exon.attributes['Name'][0],
                                    'exon_attribute':str(exon.attributes),
                                    'exon_start':exon.start,
                                    'exon_end':exon.end,
                                    'intronID':intron.attributes['ID'][0],
                                    'intron_start':intron.start,
                                    'intron_end':intron.end,
                                    'intron_attribute':str(intron.attributes),
                                    'flag_intron_retention':1},ignore_index=True)
        if hasOverlap == False:
            exonDF = exonDF.append({'chr':exon.chrom,
                                    'exonID':exon.attributes['Name'][0],
                                    'exon_attribute':str(exon.attributes),
                                    'exon_start':exon.start,
                                    'exon_end':exon.end,
                                    'intronID':np.nan,
                                    'intron_start':np.nan,
                                    'intron_end':np.nan,
                                    'intron_attribute':np.nan,
                                    'flag_intron_retention':0},ignore_index=True)
    
    # Get unique exons, flags, and sets of intron IDs they overlap with
    exonDF = exonDF.fillna(-1)
    exonDF2 = exonDF.groupby('exonID').agg({'exon_attribute':['first'],
                                            'intronID':[lambda x: "|".join(x.map(str))],
                                            'intron_attribute':['max'],
                                            'flag_intron_retention':['max']}).reset_index()
    exonDF2.columns = exonDF2.columns.droplevel(1)
    exonDF2['intronID'] = np.where(exonDF2['intron_attribute']==-1,np.nan,exonDF2['intronID'])
    del(exonDF2['intron_attribute'])
    print("{} out of {} exons have intron retention event".format(exonDF2['flag_intron_retention'].sum(),len(exonDF2)))
    
    # Merge flags with exon fragmnet IDs
    fragDF = pd.read_csv(args.inFrag,low_memory=False)
    splitDF = split_column_by_sep(df=fragDF,col_name='exon_id',sep="|")
    del(fragDF)
    mergeDF = pd.merge(splitDF,exonDF2,how='outer',left_on='exon_id',
                       right_on='exonID',indicator='merge_check')
    if mergeDF['merge_check'].value_counts()['right_only'] != 0 or mergeDF['merge_check'].value_counts()['left_only'] != 0:
        print('ERROR: Improper merge of exon fragments and intron retention flags')
    del(splitDF)
    
    # Get unique fragments
    exonNum = mergeDF['exon_id'].str.split(':').str[1].astype(int)
    sortedMergeDF = mergeDF.iloc[(mergeDF['exon_id'].map({}).fillna(exonNum)).argsort()]
    del(exonNum)
    mergeDF2 = sortedMergeDF.groupby('fragment_id').agg(
            {'chr':['first'],
             'fragment_start':['first'],
             'fragment_stop':['first'],
             'fusion_id':['first'],
             'fusion_start':['first'],
             'fusion_stop':['first'],
             'num_exons_per_frag':['first'],
             'transcript_id':['first'],
             'num_xscripts_per_frag':['first'],
             'gene_id':['first'],
             'num_genes_per_frag':['first'],
             'annotation_frequency':['first'],
             'flag_multigene':['first'],
             'total_transcripts_per_genes':['first'],
             'exon_id':[lambda x: "|".join(x.map(str))],
             'intronID':[lambda x: "|".join(x.map(str))],
             'flag_intron_retention':'max'}).reset_index()
    mergeDF2.columns = mergeDF2.columns.droplevel(1)
    del(sortedMergeDF)
    del(mergeDF)
    print("{} out of {} fragments are associated with an intron retention exon".format(
            mergeDF2['flag_intron_retention'].sum(),len(mergeDF2)))

    # Output file intron retention flagged fragments
    mergeDF2.to_csv("{}/{}_exon_fragment_annotations_flag_IR.csv".format(
            args.outDir,args.outPrefix),index=False)

if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
