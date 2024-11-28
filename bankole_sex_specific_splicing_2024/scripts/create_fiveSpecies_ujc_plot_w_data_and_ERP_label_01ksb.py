#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 12:52:24 2024

@author: k.bankole

DEVELOPED BASED ON CODE FROM APY (THANK YOU APY!!!)
"""

import pandas as pd
import numpy as np

import argparse

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import ScalarFormatter
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

from upsetplot import UpSet


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Plots the fiveSpecies UJC + other info and labels. Input: "
        "fiveSpecies UJC GTF (--fs-file), "

        "self mapped ujc_xscript_link (--sm-link), "

        "datafile with raw cnts of reads per UJC per sample ( --data-file), "

        "datacol (--data-col), which is a list of the columns which "
        "will be used to calculate prop_F (these cols will "
        "be averaged across sample to create a number of reads "
        "for male/female -> prop_F), "

        "annotated infoERP file (--info-erp-file), "

        "ER GTF (--er-file) and ES GTF (--es-file), "

        "fiveSpecies flag file (--flag-file), "

        "desired input gene to subset input (--gene), "

        "genome of files (--genome). "

        "Will output 2 plots: plot of UJCs per gene with: ERP labelled, "
        "gene model exon regions/exon segments, UJCs colored based on the "
        "proportion of reads associated with the jxnHash that are female vs male. "
        "Also outputs an upset plot for that gene that shows the annotations that "
        "contributed UJCs for the fiveSpecies on this genome."
    )

    # Input data
    parser.add_argument(
        "-fs",
        "--fs-file",
        dest="fsFile",
        required=True,
        help="Path to fiveSpecies_2_genome_ujc.gtf"
    )

    parser.add_argument(
        "-sm",
        "--sm-link",
        dest="smLinkFile",
        required=True,
        help="Path to self-mapped link file "
        "(annoName_2_genome_ujc_xscript_link.csv)"
    )

    parser.add_argument(
        "-d",
        "--data-file",
        dest="dataFile",
        required=True,
        help="Path to datafile_jxnHash_genome_anno.csv"
    )

    parser.add_argument(
        "-c",
        "--data-col",
        nargs="+",
        dest="dataCol",
        required=True,
        help="List the dataCol these space separated."
    )

    parser.add_argument(
        "-erp",
        "--info-erp-file",
        dest="erpFile",
        required=True,
        help="Path to fiveSpecies_2_genome_ujc_er_vs_fiveSpecies_2_genome_ujc_infoERP.csv"
    )

    parser.add_argument(
        "-er",
        "--er-file",
        dest="erFile",
        required=True,
        help="Path to fiveSpecies_2_genome_ujc_er.gtf"
    )

    parser.add_argument(
        "-es",
        "--es-file",
        dest="esFile",
        required=True,
        help="Path to fiveSpecies_2_genome_ujc_es.gtf"
    )

    parser.add_argument(
        "-f",
        "--flag-file",
        dest="flagFile",
        required=True,
        help="Path to flag_fiveSpecies_2_genome_ujc.csv"
    )

    parser.add_argument(
        "-gn",
        "--gene",
        dest="gene",
        required=True,
        help="Desired gene"
    )

    parser.add_argument(
        "-sym",
        "--symbol",
        dest="geneSymbol",
        required=False,
        help="Gene symbol (not required)"
    )

    parser.add_argument(
        "-gm",
        "--genome",
        dest="genome",
        required=True,
        help="Genome"
    )

    # Output data
    parser.add_argument(
        "-o",
        "--output-directory",
        dest="outdir",
        required=True,
        help="Path to output directory for plots (2 files)"
    )

    args = parser.parse_args()
    return args


def readExonData(gtfFile):

    # Read in a GTF. Create geneID and transcriptID column based on attributes.

    columns = ['seqname', 'source', 'feature', 'start',
               'end', 'score', 'strand', 'frame', 'attribute']

    data = pd.read_csv(gtfFile, sep='\t', comment='#',
                       header=None, names=columns,
                       low_memory=False)

    data = data[data['feature'] == "exon"]

    data['start'] = data['start'].astype(int)
    data['end'] = data['end'].astype(int)
    data['geneID'] = data['attribute'].str.extract(r'gene_id "([^"]+)"')
    data['transcriptID'] = data['attribute'].str.extract(
        r'transcript_id "([^"]+)"')

    data.drop(columns=['score', 'frame', 'attribute'], inplace=True)

    return data


def createFSPanel(fsDfr, plotInfoDfr, fsPanel, exonHeight, alpha, colorLst, longestTrLngth):
    # # Used for testing the colors. Makes a range of props from 0.1 to 1
    # fsPlottingDfr['prop_F'] = [0.1 * i for i in range(len(fsPlottingDfr))]

    # Scale padding based on the number of UJC in gene and length of UJCs
    xPadding = 0.1 * longestTrLngth
    yPadding = .005 * len(plotInfoDfr)

    # Loop through every jxnhash in the fiveSpecies (every row)
    for row in plotInfoDfr.to_dict('records'):

        # Extract variables from row in dataframe
        figurePos = row['figurePos']
        exonLst = row['coords']

        dataProp = row['prop_F']

        isSelfMap = row['isSelfMap']

        # Set up variables for coloring UJCs based on proportion of reads that are F
        colorMap = LinearSegmentedColormap.from_list(
            "female_v_male_cmap", colorLst)
        normalize = Normalize(vmin=0, vmax=1)

        # Loop through every exon for that jxnHash in the fiveSpecies
        prevEnd = None
        for exon in exonLst:

            start = exon[0]
            end = exon[1]

            # Set the color of the rectangles and lines.
            # Black = 0 reads (not in data)
            # Red (100% female) to Purple (50/50) to Blue (100% male)
            if dataProp == np.nan:
                color = 'black'
            else:
                color = colorMap(normalize(dataProp))

            # Add a rectangle to figure for every exon using start/end pos
            fsPanel.add_patch(Rectangle(
                xy=(start, figurePos - exonHeight / 2),
                width=end - start,
                height=exonHeight,
                facecolor=color,
                edgecolor=color,
                alpha=alpha
            ))

            # After the first exon, plot a line from the start of the last exon to the beginning of the first
            if prevEnd is not None:
                fsPanel.plot(
                    [prevEnd, start],
                    [figurePos, figurePos],
                    color=color,
                    alpha=alpha
                )

            prevEnd = end

        # If the UJC is in the self-mapped as well, create a small horizontal line next to transcript
        # Can be changed from horizontal line later

        if isSelfMap:
            fsPanel.hlines(
                y=figurePos,
                xmin=fsDfr['start'].min() - (.7 * xPadding),
                xmax=fsDfr['start'].min() - (.3 * xPadding),
                color='black',
                alpha=alpha
            )

    # Create an object for legend for the range of colors (will be added to plot later everything is added)
    colorLegend = ScalarMappable(norm=normalize, cmap=colorMap)
    colorLegend.set_array([])

    # # TODO: Reverse plot based on strand? Or no?
    # If yes: just swap the  xlims and use prevStart rather than prevEnd above

    # y-axis
    # Set y-axis limits to go slightly past the number of jxnHash in the fiveSpecies (descending order)

    # Goes from
    fsPanel.set_ylim(len(plotInfoDfr) + yPadding + 1, -yPadding)

    # Make sure transcripts are set to their assigned y-position (=figurePos, which is based on sorted ERP)
    fsPanel.set_yticks(plotInfoDfr['figurePos'])

    # Use ERP as the label for y-axis
    fsPanel.set_yticklabels(plotInfoDfr['ERP'])

    # Set y-axis label
    fsPanel.set_ylabel("ERP")

    # x-axis
    # Set x-axis limits to go slightly past the first start/last end for exon coords in the fiveSpecies
    fsPanel.set_xlim(fsDfr['start'].min() - xPadding,
                     fsDfr['end'].max() + xPadding)

    # Assure that scientific notation is not used and that the x-axis can fit labels
    fsPanel.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    fsPanel.ticklabel_format(style="plain", axis='x')

    # Set x-axis label
    fsPanel.set_xlabel(
        "Genomic Position")

    return colorLegend


def createGeneModelPanel(erCoordLst, esCoordLst, geneModelPanel, exonHeight, alpha, color, segmentColor, fsPanel):

    prevEnd = None
    # Loop through every ER
    for er in erCoordLst:

        start = er[0]
        end = er[1]

        # Add a rectangle to figure for every exon region using start/end pos
        geneModelPanel.add_patch(Rectangle(
            xy=(start, 1 - exonHeight / 2),
            width=end - start,
            height=exonHeight,
            facecolor=color,
            edgecolor=color,
            alpha=alpha
        ))

        # After the first exon, plot a line from the start of the last exon to the beginning of the first
        if prevEnd is not None:
            geneModelPanel.plot([prevEnd, start], [1, 1],
                                color=color, alpha=alpha)
        prevEnd = end

    # Loop through every ES and plot a vertical line at the end of every segment
    # (except the last, which is just the end of the gene)
    for es in esCoordLst[:-1]:
        end = es[1]
        geneModelPanel.vlines(
            x=end,
            ymin=1 - exonHeight / 2,
            ymax=1 + exonHeight / 2,
            color=segmentColor
        )

    # TODO: Reverse based on strand?
    # If so: just swap the  xlims and use prevStart rather than prevEnd

    # y-axis
    # Set y-axis limits to large enough to fit gene model
    geneModelPanel.set_ylim(1 - exonHeight - .3, 1 + exonHeight + .3)

    # Make y-axis blank
    geneModelPanel.axes.get_yaxis().set_ticks([])

    # x-axis
    # Set x-axis limits to be same as fiveSpecies panel
    geneModelPanel.set_xlim(fsPanel.get_xlim())

    # Assure that scientific notation is not used and that the x-axis can fit labels
    geneModelPanel.xaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    geneModelPanel.ticklabel_format(style="plain", axis='x')

    # Set x-axis label
    geneModelPanel.set_xlabel("Genomic Position")

    # Title this panel of the figure
    geneModelPanel.title.set_text("Gene Model Exon Region and Segments")


def createGeneUpsetPlot(flagDfr):

    # Transform flag file to fit into upset plot requirements
    flagDfr['jxnHash'] = flagDfr[[
        col for col in flagDfr.columns if '_jxnHash' in col]]

    flagDfr = flagDfr.rename(columns={
        'flag_dmel650_2_dmel6_ujc': 'dmel650',
        'flag_dsimWXD_2_dsim2_ujc': 'dsimWXD',
        'flag_dsim202_2_dsim2_ujc': 'dsim202',
        'flag_dsan11_2_dsan1_ujc': 'dsan11',
        'flag_dyak21_2_dyak2_ujc': 'dyak21',
        'flag_dser11_2_dser1_ujc': 'dser11',
    })

    flagDfr = flagDfr[['dser11', 'dyak21', 'dsan11',
                       'dsim202', 'dsimWXD', 'dmel650', 'jxnHash']]

    flagDfr.replace({1: True, 0: False}, inplace=True)

    flagDfr.set_index(
        [col for col in flagDfr.columns if 'jxnHash' not in col], inplace=True)

    # Create upset object
    upset = UpSet(
        flagDfr,
        subset_size="count",
        show_counts=True,
        sort_by="degree",
        sort_categories_by=None,
    )

    return upset


def main():
    fsFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc.gtf"
    smLinkFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/annotation_id_ujc_output/dmel650_2_dmel6_ujc_xscript_link.csv"
    dataFile = "/nfshome/k.bankole/mclab/SHARE/McIntyre_Lab/sex_specific_splicing/datafile_jxnHash_dmel6_anno.csv"
    dataCol = ['rawCnts_dmel_F_rep4', 'rawCnts_dmel_F_rep5', 'rawCnts_dmel_F_rep6',
               'rawCnts_dmel_M_rep4', 'rawCnts_dmel_M_rep5', 'rawCnts_dmel_M_rep6']
    erpFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er_vs_fiveSpecies_2_dmel6_ujc_infoERP.csv"
    erFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_er.gtf"
    esFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/fiveSpecies_2_dmel6_ujc_es.gtf"
    flagFile = "/nfshome/k.bankole/mnt/exasmb.rc.ufl.edu-blue/mcintyre/share/sex_specific_splicing/fiveSpecies_annotations/flag_fiveSpecies_2_dmel6_ujc.csv"
    inputGn = "FBgn0264270"
    genome = "dmel6"
    outdir = "/nfshome/k.bankole/Desktop/test_folder"
    symbol = "Sxl"
    # symbol = None

    fsFile = args.fsFile
    smLinkFile = args.smLinkFile
    dataFile = args.dataFile
    dataCol = args.dataCol
    erpFile = args.erpFile
    erFile = args.erFile
    esFile = args.esFile
    flagFile = args.flagFile
    inputGn = args.gene
    genome = args.genome
    outdir = args.outdir
    symbol = args.geneSymbol

    # Read in all input files
    inFSDfr = readExonData(fsFile)
    inSMLinkDfr = pd.read_csv(smLinkFile, low_memory=False)
    inDataDfr = pd.read_csv(dataFile, low_memory=False)
    inERPDfr = pd.read_csv(erpFile, low_memory=False)
    inERDfr = readExonData(erFile)
    inESDfr = readExonData(esFile)
    inFlagDfr = pd.read_csv(flagFile, low_memory=False)

    # Subset every input to the desired input gene
    # TODO: check that these aren't empty
    fsDfr = inFSDfr[inFSDfr['geneID'] == inputGn].copy().rename(
        columns={'transcriptID': 'jxnHash'})

    print(
        f"Number of fiveSpecies UJC in {inputGn}:", fsDfr['jxnHash'].nunique())

    longestTrLngth = fsDfr['end'].max() - fsDfr['start'].min()

    smLinkDfr = inSMLinkDfr[inSMLinkDfr['geneID'] == inputGn].copy()[
        ['transcriptID', 'jxnHash']]

    dataDfr = inDataDfr[inDataDfr['geneID'] == inputGn].copy()

    erDfr = inERDfr[inERDfr['geneID'] == inputGn].copy()
    esDfr = inESDfr[inESDfr['geneID'] == inputGn].copy()

    erpDfr = inERPDfr[inERPDfr['geneID'] == inputGn].copy()[
        ['jxnHash', 'ERP_plus', 'ERP']]

    flagDfr = inFlagDfr[inFlagDfr['geneID'] == inputGn].copy()

    # Grab strand
    strand = fsDfr['strand'].iloc[0]

    # Create DataFrame for plotting
    # DataFrame containing:
    # fiveSpecies jxnHash and exon coords
    # isSelfMap,
    # ERP (mergeERPDfr)
    # prop_F from the data used to color the fiveSpecies jxnHash

    # make fiveSpecies unique on jxnHash based on GTF with exons stored
    fsDfr['coords'] = fsDfr.apply(
        lambda row: (row['start'], row['end']), axis=1)
    fsUniqJxnHashDfr = fsDfr.groupby(['geneID', 'jxnHash']).agg(
        list)['coords'].reset_index().sort_values(by='jxnHash')

    # sort by first start
    fsUniqJxnHashDfr['coords'] = fsUniqJxnHashDfr['coords'].apply(
        lambda x: sorted(x, key=lambda y: y[1]))

    # Merge in self-mapped transcripts
    mergeSMDfr = pd.merge(fsUniqJxnHashDfr, smLinkDfr, on='jxnHash',
                          how='outer', indicator='merge_check')

    if (mergeSMDfr['merge_check'] != 'right_only').all():
        mergeSMDfr['isSelfMap'] = mergeSMDfr['merge_check'].apply(
            lambda x: x == 'both')
        mergeSMDfr.drop(['transcriptID', 'merge_check'], axis=1, inplace=True)
    else:
        raise Exception("An error occurred when merging the xscript link and fiveSpecies GTF. "
                        "There are jxnHash that only in the fiveSpecies (not possible).")

    # Get desired data read count columns to take an average of
    colorMapDfr = dataDfr[['jxnHash'] + dataCol].copy()

    # Get average number of F read across desired samples
    colorMapDfr['total_F_read'] = colorMapDfr[colorMapDfr.filter(
        like='_F_rep').columns].sum(axis=1)
    colorMapDfr['avg_F_read'] = colorMapDfr['total_F_read'] / \
        len(colorMapDfr.filter(like='_F_rep').columns)

    # Get average number of M read across desired samples
    colorMapDfr['total_M_read'] = colorMapDfr[colorMapDfr.filter(
        like='_M_rep').columns].sum(axis=1)
    colorMapDfr['avg_M_read'] = colorMapDfr['total_M_read'] / \
        len(colorMapDfr.filter(like='_M_rep').columns)

    # Create a proportion of female reads using above averages
    colorMapDfr['sumAvgRead'] = colorMapDfr['avg_F_read'] + \
        colorMapDfr['avg_M_read']

    colorMapDfr['prop_F'] = colorMapDfr.apply(
        lambda row: row['avg_F_read'] / row['sumAvgRead'] if row['sumAvgRead'] > 0 else None, axis=1)

    colorMapDfr = colorMapDfr[['jxnHash', 'prop_F']]

    # Merge in data prop_F column
    mergeDataDfr = pd.merge(mergeSMDfr, colorMapDfr,
                            on='jxnHash', how='outer', indicator='merge_check')

    # Verify there are no jxnHash that are only in the data
    # its ok for there to be anno only jxnHash (just means 0 reads for it)
    if (mergeDataDfr['merge_check'] != 'right_only').all():
        mergeDataDfr.drop('merge_check', axis=1, inplace=True)
    else:
        raise Exception("An error occurred when merging the datafile and fiveSpecies GTF. "
                        "There are jxnHash that are only in the datafile.")

    # Merge in ERPs
    mergeERPDfr = pd.merge(mergeDataDfr, erpDfr, on='jxnHash',
                           how='outer', indicator='merge_check')

    # Verify that all jxnHash are in both files
    if (mergeERPDfr['merge_check'] == 'both').all():
        mergeERPDfr.drop('merge_check', axis=1, inplace=True)
    else:
        raise Exception("An error occurred when merging the ERP and fiveSpecies GTF. "
                        "There are jxnHash that aren't in both.")

    # Sort by ERP (this dfr is used for plotting FS)
    plotInfoDfr = mergeERPDfr.sort_values(
        by='ERP', ignore_index=True)

    # Create a figure y-position for jxnHash based on sorted ERP
    plotInfoDfr = plotInfoDfr.reset_index(names='figurePos')
    plotInfoDfr['figurePos'] = plotInfoDfr['figurePos'] + 1

    # 3. Create lists of the ERs and ESs for the gene

    # make er and es gtf unique on geneID based on GTF with exons stored
    erDfr['ER_coords'] = erDfr.apply(
        lambda row: (row['start'], row['end']), axis=1)
    geneERDfr = erDfr.groupby(['geneID']).agg(
        list)['ER_coords'].reset_index()

    esDfr['ES_coords'] = esDfr.apply(
        lambda row: (row['start'], row['end']), axis=1)
    geneESDfr = esDfr.groupby(['geneID']).agg(
        list)['ES_coords'].reset_index()

    # TODO: check merge??? doesn't seem necessary here since its just one gene
    # Create dict with ER_coords and ES_coords
    geneERESCoordDct = pd.merge(geneERDfr, geneESDfr, on='geneID',
                                how='outer').set_index('geneID').iloc[0].to_dict()

    # Use dict to create a list of ERs and ESs (use for plotting gene model)
    erCoordLst = geneERESCoordDct['ER_coords']
    esCoordLst = geneERESCoordDct['ES_coords']

    # Setup Plot and Panels
    fig = plt.figure(figsize=(15, 10))
    fsPanel = plt.subplot(1, 1, 1)

    # Plot FS UJC and return color legend for creating colorbar legend later
    colorLegend = (
        createFSPanel(
            fsDfr=fsDfr,
            plotInfoDfr=plotInfoDfr,
            fsPanel=fsPanel,
            exonHeight=0.6,
            alpha=0.6,
            colorLst=['blue', 'purple', 'red'],
            longestTrLngth=longestTrLngth
        ))

    # Plot gene model underneath fs panel
    divider = make_axes_locatable(fsPanel)
    geneModelPanel = divider.append_axes('bottom', size="20%")

    createGeneModelPanel(
        erCoordLst=erCoordLst,
        esCoordLst=esCoordLst,
        geneModelPanel=geneModelPanel,
        exonHeight=0.6,
        alpha=0.6,
        color='lightblue',
        segmentColor='black',
        fsPanel=fsPanel,
    )

    fig.tight_layout()

    # Create separate panel for color bar and move it into a good position
    # left, bottom, width, height
    cbPanel = fig.add_axes([.6, 0.353, 0.5, 0.515], frameon=False)
    cbPanel.set_xticks([])
    cbPanel.set_yticks([])
    colorBar = plt.colorbar(colorLegend, ax=cbPanel, orientation="vertical")
    labelText = ("Proportion of Female Reads in Data\n"
                 "(0 = 100% M, 1 = 100% F,\n "
                 "Black = Not Supported by Data)")
    colorBar.set_label(labelText)
    # Center align the label text
    colorBar.ax.xaxis.label.set_horizontalalignment('center')

    # Title the whole figure

    if symbol:
        fsPanel.title.set_text(f"{genome} - {inputGn} - {symbol} (strand={strand})\n"
                               "fiveSpecies Annotated UJC\n\n"
                               f"Number of UJCs: {len(plotInfoDfr)}")
    else:
        fsPanel.title.set_text(f"{genome} - {inputGn} (strand={strand})\n"
                               "fiveSpecies Annotated UJC\n\n"
                               f"Number of UJCs: {len(plotInfoDfr)}")

    # Create the upset plot (still need an output step)
    upset = createGeneUpsetPlot(flagDfr)
    fig2 = plt.figure(figsize=(23, 7))
    upset.plot(fig=fig2)['totals'].set_title("# UJCs in OG Anno")

    fig2.tight_layout()
    plt.ylabel('Number of UJCs in OG Anno')

    if symbol:
        fig2.suptitle(f"{genome} - {inputGn} - {symbol}:\n"
                      "Source of fiveSpecies UJCs when aligned to {genome}")
    else:
        fig2.suptitle(f"{genome} - {inputGn}:\n"
                      "Source of fiveSpecies UJCs when aligned to {genome}")

    # Save both plots

    if symbol:
        plotFile = f"{outdir}/plot_combo_ujc_data_erp_{inputGn}_{symbol}_fiveSpecies_2_{genome}.svg"
        upsetFile = f"{outdir}/upst_ujc_anno_source_{inputGn}_{symbol}_fiveSpecies_2_{genome}.svg"
    else:
        plotFile = f"{outdir}/plot_combo_ujc_data_erp_{inputGn}_fiveSpecies_2_{genome}.svg"
        upsetFile = f"{outdir}/upst_ujc_anno_source_{inputGn}_fiveSpecies_2_{genome}.svg"

    fig.savefig(plotFile, dpi=600, format="svg", bbox_inches='tight')
    print(f"Saved plot: {plotFile}")

    fig2.savefig(upsetFile, dpi=600, format="svg", bbox_inches='tight')
    print(f"Saved plot: {upsetFile}")

    # fig.show()
    # fig2.show()

    plt.close(fig)
    plt.close(fig2)


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
