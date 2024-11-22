#!/usr/bin/env python

import pandas as pd
import time
import os.path


# DESCRIPTION: Create gene model info table for each species. Unique on gene. Columns:
# num transcripts per gene in the self-mapped GTF (or both in the case of dsim2)
# num UJC per gene in fiveSpecies
# num unique ERP per gene when doing fiveSpecies ER vs UJC
# num exon region per gene
# num exon segment per gene

# Stolen from TranD io
def read_exon_data_from_file(infile, keepSrc=False):
    """
    Create a pandas dataframe with exon records from a gtf file
    Raw gtf data:
    seqname source  feature  start end  score strand frame attributes
    2L FlyBase 5UTR  7529  7679 . +  . gene_symbol "CG11023"; transcript_id "FBtr...
    2L FlyBase exon  7529  8116 . +  . gene_symbol "CG11023"; transcript_id "FBtr...

    Exon Fragment Data for Event Analysis:
    seqname start end   strand gene_id      transcript_id
    2L      7529  8116  +      FBgn0031208  FBtr0300689
    2L      8193  9484  +      FBgn0031208  FBtr0300689
    2L      9839  11344 -      FBgn0002121  FBtr0078169
    """

    print("Reading GTF...")

    omegatic = time.perf_counter()

    all_gtf_columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame',
                       'attributes', 'comments']

    drop_columns = ['feature', 'score', 'frame', 'comments']
    # drop_columns = ['source', 'feature', 'score', 'frame', 'comments']

    data = pd.read_csv(infile, sep='\t', comment='#',
                       header=None, low_memory=False)
    file_cols = data.columns

    if len(file_cols) < len(all_gtf_columns):
        gtf_cols = all_gtf_columns[:len(file_cols)]
    data.columns = gtf_cols
    drop_cols = [x for x in drop_columns if x in gtf_cols]

    data = data[data['feature'] == 'exon']
    data = data.drop(labels=drop_cols, axis=1)

    data['source'] = data['source'].astype(str)
    data['seqname'] = data['seqname'].astype(str)
    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)

    data.reset_index(drop=True, inplace=True)

    sourceLst = []
    seqnameLst = []
    startLst = []
    endLst = []
    strandLst = []
    geneIDLst = []
    xscriptIDLst = []

    for row in data.to_dict('records'):
        rawAttr = row['attributes']
        attrLst = [x.strip() for x in rawAttr.strip().split(';')]
        gnTrAttr = [
            x for x in attrLst if 'transcript_id' in x or 'gene_id' in x]
        gene_id, transcript_id = None, None

        for item in gnTrAttr:
            if 'gene_id' in item:
                gene_id = item.split('gene_id')[1].strip().strip('\"')
            elif 'transcript_id' in item:
                transcript_id = item.split(
                    'transcript_id')[-1].strip().strip('\"')

        if not gene_id:
            print("gene_id not found in '{}'", row)
            gene_id = None

        if not transcript_id:
            print("transcript_id not found in '{}'", row)
            transcript_id = None

        sourceLst.append(row['source'])
        seqnameLst.append(row['seqname'])
        startLst.append(row['start'])
        endLst.append(row['end'])
        strandLst.append(row['strand'])

        geneIDLst.append(gene_id)
        xscriptIDLst.append(transcript_id)

    if keepSrc:
        newData = pd.DataFrame(
            {
                'seqname': seqnameLst,
                'source': sourceLst,
                'start': startLst,
                'end': endLst,
                'strand': strandLst,
                'gene_id': geneIDLst,
                'transcript_id': xscriptIDLst
            })
    else:

        newData = pd.DataFrame(
            {
                'seqname': seqnameLst,
                'start': startLst,
                'end': endLst,
                'strand': strandLst,
                'gene_id': geneIDLst,
                'transcript_id': xscriptIDLst
            })

    print("Exon data rows:", newData.shape[0])

    missing_value_num = newData.isnull().sum().sum()
    if missing_value_num > 0:
        print("Total number of missing values:", missing_value_num)
    else:
        print("No missing values in data")

    gene_id_missing_value_num = newData['gene_id'].isnull().sum()

    transcript_id_missing_value_num = newData['transcript_id'].isnull().sum()

    if gene_id_missing_value_num > 0:
        print("Missing gene_id value number:",
              gene_id_missing_value_num)
    if transcript_id_missing_value_num > 0:
        print("Missing transcript_id value number:",
              transcript_id_missing_value_num)

    newData['start'] = pd.to_numeric(newData['start'], downcast="unsigned")
    newData['end'] = pd.to_numeric(newData['end'], downcast="unsigned")

    toc = time.perf_counter()

    print(f"GTF Read complete,  took {toc-omegatic:0.4f} seconds.")
    return newData


# Change if necessary
PROJ = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing"

# Genomes and their self-mapped GTF
genomeRefDct = {

    'dmel6': ['dmel_fb650/dmel650_2_dmel6.gtf'],
    'dsim2': ['dsim_fb202/dsim202_2_dsim2.gtf', 'dsim_fb202/dsimWXD_2_dsim2.gtf'],
    'dsan1': ['dsan_Prin_1.1/dsan11_2_dsan1.gtf'],
    'dyak2': ['dyak_Prin_Tai18E2_2.1/dyak21_2_dyak2.gtf'],
    'dser1': ['dser1.1/dser11_2_dser1.gtf']
}


# SM == SELF-MAPPED
# FS == FIVESPECIES

# TODO: add geneSymbol/mel_geneSymbol??
for genome, referenceLst in genomeRefDct.items():

    print(f"Creating table for {genome}...")

    # Define input files
    inSMGTFLst = [
        '/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/references/' + x for x in referenceLst]
    inFSGTF = f"{PROJ}/fiveSpecies_annotations/fiveSpecies_2_{genome}_ujc.gtf"
    inERPFile = f"{PROJ}/fiveSpecies_annotations/fiveSpecies_2_{genome}_ujc_er_vs_fiveSpecies_2_{genome}_ujc_infoERP.csv"
    inERGTF = f"{PROJ}/fiveSpecies_annotations/fiveSpecies_2_{genome}_ujc_er.gtf"
    inESGTF = f"{PROJ}/fiveSpecies_annotations/fiveSpecies_2_{genome}_ujc_es.gtf"

    # Set up empty variable for Dfr
    smTrPerGeneDfr = pd.DataFrame()

    # Loop through self-mapped gtfs (done because dsim2 has 2)
    for gtfFile in inSMGTFLst:

        # Get name of GTF for column name later
        gtfName = os.path.basename(gtfFile).split('.gtf')[0]

        # Read in sm GTF
        gtfDfr = read_exon_data_from_file(gtfFile)

        # If there is only one sm GTF (dfr is empty),
        # -> create new Dfr that are unique on gene w num sm transcript gene
        if smTrPerGeneDfr.empty:
            smTrPerGeneDfr = gtfDfr.groupby(
                'gene_id')['transcript_id'].nunique().reset_index()
            smTrPerGeneDfr.columns = ['geneID', f'num_{gtfName}_transcript']
        else:
            # if there is more than one sm GTF (dfr is not empty, loop occurs more than once)
            # merge in the sm transcript per gene from the additional sm gtf(s)

            tempTrPerGeneDfr = gtfDfr.groupby(
                'gene_id')['transcript_id'].nunique().reset_index()
            tempTrPerGeneDfr.columns = ['geneID', f'num_{gtfName}_transcript']

            # There is going to be left and right only because not all annotations cover genes the same
            # (ex: some dsim2 genes have no transcripts in dsim202 after mapping but have some in dsimWXD)
            # -> no merge_check
            smTrPerGeneDfr = pd.merge(
                smTrPerGeneDfr, tempTrPerGeneDfr, on='geneID', how='outer').fillna(0)

            del tempTrPerGeneDfr

    # Read in FS GTF
    fsGTFDfr = read_exon_data_from_file(inFSGTF)

    # -> create new Dfr that is unique on gene w num ujc per gene

    fsTrPerGeneDfr = fsGTFDfr.groupby(
        'gene_id')['transcript_id'].nunique().reset_index()
    fsTrPerGeneDfr = fsTrPerGeneDfr.rename(
        columns={'gene_id': 'geneID', 'transcript_id': 'num_fiveSpecies_UJC'})

    # Read in ERP info file for fiveSpecies ER vs fiveSpecies UJC
    erpDfr = pd.read_csv(inERPFile, low_memory=False)

    # -> create new Dfr that is unique on gene w num unique ERP per gene
    erpPerGeneDfr = erpDfr.groupby(
        'geneID')['ERP_plus'].nunique().reset_index()
    erpPerGeneDfr = erpPerGeneDfr.rename(
        columns={'ERP_plus': 'num_unique_ERP'})

    # Read in ER and ES GTF
    erGTFDfr = read_exon_data_from_file(inERGTF)
    esGTFDfr = read_exon_data_from_file(inESGTF)

    # -> create new Dfrs that are unique on gene w num ER and ES per gene
    # simply counts number of exon rows per gene (verified that this correctly counts num ER/ES)
    erPerGeneDfr = erGTFDfr.groupby('gene_id').size().reset_index()
    erPerGeneDfr = erPerGeneDfr.rename(
        columns={'gene_id': 'geneID', 0: 'num_exonRegion'})
    esPerGeneDfr = esGTFDfr.groupby('gene_id').size().reset_index()
    esPerGeneDfr = esPerGeneDfr.rename(
        columns={'gene_id': 'geneID', 0: 'num_exonSegment'})

    # MERGING!!!!

    # MERGE FS TRANSCRIPTS PER GENE AND SM TRANSCRIPTS PER GENE
    # It is ok for there to be left_only and right_only here -> no merge_check
    # left_only: most likely a gene that had all of its transcripts switched to a different gene due to GFFCompare
    # right_only: genes that have no self-mapped transcripts, but have UJCs from other species when mapped to mel
    mergeDf = pd.merge(smTrPerGeneDfr, fsTrPerGeneDfr,
                       on='geneID', how='outer').fillna(0)

    # MERGE IN NUM UNIQ ERP PER GENE
    mergeDf = pd.merge(mergeDf, erpPerGeneDfr, on='geneID',
                       how='outer', indicator='merge_check')

    # Check for no right_only
    if (mergeDf['merge_check'] == 'right_only').any():
        raise Exception(
            "There was an issue when merging in the number of ERPs per gene. There are genes that are only in the ERP output.")
    else:
        mergeDf.drop('merge_check', axis=1, inplace=True)

    # MERGE IN NUM ER
    mergeDf = pd.merge(mergeDf, erPerGeneDfr, on='geneID',
                       how='outer', indicator='merge_check')

    # Check for no right_only
    if (mergeDf['merge_check'] == 'right_only').any():
        raise Exception(
            "There was an issue when merging in the number of ERPs per gene. There are genes that are only in the ERP output.")
    else:
        mergeDf.drop('merge_check', axis=1, inplace=True)

    # MERGE IN NUM ES
    mergeDf = pd.merge(mergeDf, esPerGeneDfr, on='geneID',
                       how='outer', indicator='merge_check')

    # Check for no right_only
    if (mergeDf['merge_check'] == 'right_only').any():
        raise Exception(
            "There was an issue when merging in the number of ERPs per gene. There are genes that are only in the ERP output.")
    else:
        mergeDf.drop('merge_check', axis=1, inplace=True)

    mergeDf = mergeDf.fillna(0)
    outfile = f"/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/table_{genome}_gene_model_info.csv"
    mergeDf.to_csv(outfile, index=False)
