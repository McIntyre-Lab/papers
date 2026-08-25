#!/usr/bin/env python

from upsetplot import UpSet
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# TODO: CHANGE OUTPUT FILE NAMES


def plotSharedJxnHashPerGenome(flagFile, outdir):
    """
    Plot, for an entire fiveSpecies UJC GTF, an upset plot of the 
    num jxnHash against where they came from (aka a venn diagram 
    for every annotation used to create the fiveSpecies). Uses the 
    fiveSpecies flag file.

    Parameters
    ----------
    flagFile : String
        Input the path to a flag_fiveSpecies_2_*_ujc.csv file.
    outdir : TYPE
        Output directory for plot.
    """

    genome = flagFile.split("flag_fiveSpecies_2_")[1].split('_ujc.csv')[0]

    flagDfr = pd.read_csv(flagFile, low_memory=False)

    upsetDfr = flagDfr[[genome+'_jxnHash',] +
                       [col for col in flagDfr.columns if 'flag' in col]].copy(deep=True)

    upsetDfr = upsetDfr.rename(columns={
        'flag_dmel650_2_dmel6_ujc': 'dmel650',
        'flag_dsimWXD_2_dsim2_ujc': 'dsimWXD',
        'flag_dsim202_2_dsim2_ujc': 'dsim202',
        'flag_dsan11_2_dsan1_ujc': 'dsan11',
        'flag_dyak21_2_dyak2_ujc': 'dyak21',
        'flag_dser11_2_dser1_ujc': 'dser11',
    })

    # Need to order the columns in the opposite direction of the way you
    # want them to appear in the upset plot
    upsetDfr = upsetDfr[['dser11', 'dyak21', 'dsan11',
                         'dsim202', 'dsimWXD', 'dmel650', genome+'_jxnHash']]

    upsetDfr.replace({1: True, 0: False}, inplace=True)

    # Set index to flag columns (necessary for upset function)
    upsetDfr.set_index(
        [col for col in upsetDfr.columns if 'jxnHash' not in col], inplace=True)

    upset = UpSet(
        upsetDfr,
        subset_size="count",
        show_counts=True,
        sort_by="degree",
        sort_categories_by=None,
    )

    upset.plot()['totals'].set_title("Number of UJCs")
    plt.ylabel('Shared Number of UJCs')
    plt.suptitle("Number and Origin of UJCs Aligned to dmel6".format(genome))
    plt.savefig(
        outdir + f"/upst_num_shared_jxnHash_fiveSpecies_2_{genome}.png",
        dpi=600, format="png")


def plotGeneJxnHash(geneKey, flagFile, geneDct, genome, outdir):

    geneKeyDfr = pd.read_csv(geneKey, low_memory=False)
    geneKeyDfr = geneKeyDfr[['transcript_id', 'output_gene_id']]
    geneKeyDfr = geneKeyDfr.rename(columns={
        'output_gene_id': 'geneID',
        'transcript_id': 'jxnHash'
    })

    flagDfr = pd.read_csv(flagFile, low_memory=False)
    flagDfr['jxnHash'] = flagDfr[[
        col for col in flagDfr.columns if '_jxnHash' in col]]

    # merge_check is all "both!"
    mergeDfr = pd.merge(geneKeyDfr, flagDfr, on='jxnHash',
                        how='outer', indicator='merge_check')

    if (mergeDfr['merge_check'] != "both").any():
        raise Exception("Merge error.")
    else:
        mergeDfr = mergeDfr.drop(columns='merge_check')

    if (mergeDfr['geneID_x'] != mergeDfr['geneID_y']).any():
        raise Exception("Merge error.")
    else:
        mergeDfr['geneID'] = mergeDfr['geneID_x']
        mergeDfr.drop(['geneID_x', 'geneID_y'], axis=1, inplace=True)

    for geneName, geneID in geneDct.items():

        Dfr = mergeDfr[mergeDfr['geneID'].isin([geneID])]

        Dfr = Dfr.rename(columns={
            'flag_dmel650_2_dmel6_ujc': 'dmel650',
            'flag_dsimWXD_2_dsim2_ujc': 'dsimWXD',
            'flag_dsim202_2_dsim2_ujc': 'dsim202',
            'flag_dsan11_2_dsan1_ujc': 'dsan11',
            'flag_dyak21_2_dyak2_ujc': 'dyak21',
            'flag_dser11_2_dser1_ujc': 'dser11',
        })

        Dfr = Dfr[['dser11', 'dyak21', 'dsan11',
                   'dsim202', 'dsimWXD', 'dmel650', 'jxnHash']]

        Dfr.replace({1: True, 0: False}, inplace=True)

        Dfr.set_index(
            [col for col in Dfr.columns if 'jxnHash' not in col], inplace=True)

        upset = UpSet(
            Dfr,
            subset_size="count",
            show_counts=True,
            sort_by="degree",
            sort_categories_by=None,
        )

        upset.plot()['totals'].set_title("Number of UJCs")
        plt.ylabel('Number of UJCs')
        plt.suptitle(
            f"{geneName}: Number of UJCs when aligned to {genome}")
        plt.savefig(
            f"{outdir}/upst_num_shared_jxnHash_{geneName}_fiveSpecies_2_{genome}.png",
            dpi=600, format="png", bbox_inches='tight')


def plotSharedERPPerGenome(flagFile, erpFile, genome, outdir):

    flagDfr = pd.read_csv(flagFile, low_memory=False)
    erpDfr = pd.read_csv(erpFile, low_memory=False)

    erpDfr = erpDfr[['jxnHash', 'ERP', 'geneID']].copy().rename(
        columns={'jxnHash': f'{genome}_jxnHash'})

    mergeDfr = pd.merge(flagDfr, erpDfr,
                        on=f'{genome}_jxnHash', how='outer',
                        indicator='merge_check')

    if (mergeDfr['merge_check'] != "both").any():

        numMssngJxnHash = len(mergeDfr[mergeDfr['merge_check'] != "both"])

        print(
            f"There are {numMssngJxnHash} jxnHashes without ERPs in {genome} (due to being removed by TranD)")
        print("WARNING THESE JXNHASHES ARE NOT INCLUDED IN THE PLOT")

        workingDfr = mergeDfr[mergeDfr['merge_check'] == "both"]
        workingDfr = mergeDfr.drop('merge_check', axis=1)

    else:
        workingDfr = mergeDfr.drop('merge_check', axis=1)

    flagCol = [col for col in workingDfr.columns if 'flag' in col]
    aggregations = {
        f'{genome}_jxnHash': 'count',
    }
    for col in flagCol:
        aggregations[col] = 'max'

    uniqOnERPDfr = workingDfr.groupby('ERP').agg(aggregations).reset_index()
    flagCol = [col for col in uniqOnERPDfr.columns if 'flag' in col]
    uniqOnERPDfr[flagCol] = uniqOnERPDfr[flagCol].astype(bool)

    upsetDfr = uniqOnERPDfr.rename(columns={
        'flag_dmel650_2_dmel6_ujc': 'dmel650',
        'flag_dsimWXD_2_dsim2_ujc': 'dsimWXD',
        'flag_dsim202_2_dsim2_ujc': 'dsim202',
        'flag_dsan11_2_dsan1_ujc': 'dsan11',
        'flag_dyak21_2_dyak2_ujc': 'dyak21',
        'flag_dser11_2_dser1_ujc': 'dser11',
    }).copy()

    upsetDfr = upsetDfr[['dser11', 'dyak21', 'dsan11', 'dsim202',
                         'dsimWXD', 'dmel650', 'ERP', f'{genome}_jxnHash']]

    upset = UpSet(
        upsetDfr.set_index(['dser11', 'dyak21', 'dsan11',
                            'dsim202', 'dsimWXD', 'dmel650']),
        subset_size="count",
        show_counts=True,
        sort_by="degree",
        sort_categories_by=None,
    )

    upset.add_catplot(
        value=f'{genome}_jxnHash',
        kind="box",
        elements=4,
        showfliers=False,
        color=sns.color_palette("colorblind", 15).as_hex()[0],
    )
    upset.plot()['totals'].set_title("Number of ERPs")
    plt.ylabel('Number of UJCs')
    plt.subplots_adjust(right=1.00001)
    plt.suptitle(f"Number of ERPs within {genome} fiveSpecies annotation")
    plt.savefig(f"{outdir}/upst_num_shared_ERP_fiveSpecies_2_{genome}.png",
                dpi=600, format="png")
