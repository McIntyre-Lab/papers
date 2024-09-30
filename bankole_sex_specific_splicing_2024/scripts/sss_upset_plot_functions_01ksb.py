#!/usr/bin/env python

from upsetplot import UpSet
import pandas as pd
import matplotlib.pyplot as plt
import glob
import seaborn as sns


def plotSharedJxnHashPerGenome(flagFile, outdir):

    genome = flagFile.split("flag_fiveSpecies_2_")[1].split('_ujc.csv')[0]

    flagDf = pd.read_csv(flagFile, low_memory=False)

    upsetDf = flagDf[[genome+'_jxnHash',] +
                     [col for col in flagDf.columns if 'flag' in col]].copy(deep=True)

    upsetDf = upsetDf.rename(columns={
        'flag_dmel650_2_dmel6_ujc': 'dmel650',
        'flag_dsimWXD_2_dsim2_ujc': 'dsimWXD',
        'flag_dsim202_2_dsim2_ujc': 'dsim202',
        'flag_dsan11_2_dsan1_ujc': 'dsan11',
        'flag_dyak21_2_dyak2_ujc': 'dyak21',
        'flag_dser11_2_dser1_ujc': 'dser11',
    })

    # upsetDf = upsetDf[['dmel650','dsimWXD','dsim202','dsan11','dyak21','dser11',genome+'_jxnHash']]
    upsetDf = upsetDf[['dser11', 'dyak21', 'dsan11',
                       'dsim202', 'dsimWXD', 'dmel650', genome+'_jxnHash']]

    upsetDf.replace({1: True, 0: False}, inplace=True)

    upsetDf.set_index(
        [col for col in upsetDf.columns if 'jxnHash' not in col], inplace=True)

    upset = UpSet(
        upsetDf,
        subset_size="count",
        show_counts=True,
        sort_by="degree",
        sort_categories_by=None,
    )

    upset.plot()['totals'].set_title("Number of UJCs")
    plt.ylabel('Number of UJCs')
    plt.suptitle("Number of UJCs when aligned to {}".format(genome))
    plt.savefig(
        outdir + f"/upst_num_shared_jxnHash_fiveSpecies_2_{genome}.png",
        dpi=600, format="png")


def plotGeneJxnHash(geneKey, flagFile, geneDct, genome, outdir):

    geneKeyDf = pd.read_csv(geneKey, low_memory=False)
    geneKeyDf = geneKeyDf[['transcript_id', 'output_gene_id']]
    geneKeyDf = geneKeyDf.rename(columns={
        'output_gene_id': 'geneID',
        'transcript_id': 'jxnHash'
    })

    flagDf = pd.read_csv(flagFile, low_memory=False)
    flagDf['jxnHash'] = flagDf[[col for col in flagDf.columns if '_jxnHash' in col]]

    # merge_check is all "both!"
    mergeDf = pd.merge(geneKeyDf, flagDf, on='jxnHash',
                       how='outer', indicator='merge_check')

    if (mergeDf['merge_check'] != "both").any():
        raise Exception("Merge error.")
    else:
        mergeDf = mergeDf.drop(columns='merge_check')

    if (mergeDf['geneID_x'] != mergeDf['geneID_y']).any():
        raise Exception("Merge error.")
    else:
        mergeDf['geneID'] = mergeDf['geneID_x']
        mergeDf.drop(['geneID_x', 'geneID_y'], axis=1, inplace=True)

    for geneName, geneID in geneDct.items():

        df = mergeDf[mergeDf['geneID'].isin([geneID])]

        df = df.rename(columns={
            'flag_dmel650_2_dmel6_ujc': 'dmel650',
            'flag_dsimWXD_2_dsim2_ujc': 'dsimWXD',
            'flag_dsim202_2_dsim2_ujc': 'dsim202',
            'flag_dsan11_2_dsan1_ujc': 'dsan11',
            'flag_dyak21_2_dyak2_ujc': 'dyak21',
            'flag_dser11_2_dser1_ujc': 'dser11',
        })

        df = df[['dser11', 'dyak21', 'dsan11',
                 'dsim202', 'dsimWXD', 'dmel650', 'jxnHash']]

        df.replace({1: True, 0: False}, inplace=True)

        df.set_index(
            [col for col in df.columns if 'jxnHash' not in col], inplace=True)

        upset = UpSet(
            df,
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

    flagDf = pd.read_csv(flagFile, low_memory=False)
    erpDf = pd.read_csv(erpFile, low_memory=False)

    erpDf = erpDf[['jxnHash', 'ERP', 'geneID']].copy().rename(
        columns={'jxnHash': f'{genome}_jxnHash'})

    mergeDf = pd.merge(flagDf, erpDf,
                       on=f'{genome}_jxnHash', how='outer',
                       indicator='merge_check')

    if (mergeDf['merge_check'] != "both").any():

        numMssngJxnHash = len(mergeDf[mergeDf['merge_check'] != "both"])

        print(
            f"There are {numMssngJxnHash} jxnHashes without ERPs in {genome} (due to being removed by TranD)")
        print("WARNING THESE JXNHASHES ARE NOT INCLUDED IN THE PLOT")

        workingDf = mergeDf[mergeDf['merge_check'] == "both"]
        workingDf = mergeDf.drop('merge_check', axis=1)

    else:
        workingDf = mergeDf.drop('merge_check', axis=1)

    flagCol = [col for col in workingDf.columns if 'flag' in col]
    aggregations = {
        f'{genome}_jxnHash': 'count',
    }
    for col in flagCol:
        aggregations[col] = 'max'

    uniqOnERPDf = workingDf.groupby('ERP').agg(aggregations).reset_index()
    flagCol = [col for col in uniqOnERPDf.columns if 'flag' in col]
    uniqOnERPDf[flagCol] = uniqOnERPDf[flagCol].astype(bool)

    upsetDf = uniqOnERPDf.rename(columns={
        'flag_dmel650_2_dmel6_ujc': 'dmel650',
        'flag_dsimWXD_2_dsim2_ujc': 'dsimWXD',
        'flag_dsim202_2_dsim2_ujc': 'dsim202',
        'flag_dsan11_2_dsan1_ujc': 'dsan11',
        'flag_dyak21_2_dyak2_ujc': 'dyak21',
        'flag_dser11_2_dser1_ujc': 'dser11',
    }).copy()

    upsetDf = upsetDf[['dser11', 'dyak21', 'dsan11', 'dsim202',
                       'dsimWXD', 'dmel650', 'ERP', f'{genome}_jxnHash']]

    upset = UpSet(
        upsetDf.set_index(['dser11', 'dyak21', 'dsan11',
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
