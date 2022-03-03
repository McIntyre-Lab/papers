#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import math


def getoptions():
    parser = argparse.ArgumentParser(description='calculate and visualize the distribution of feature number statistics from a standard file input')
    parser.add_argument("-ig", "--input_g", dest="input_g", required=True, help="Input gene-transcript pairs file name")
    parser.add_argument("-ie", "--input_e", dest="input_e", required=True, help="Input exon region file name")
    parser.add_argument("-if", "--input_f", dest="input_f", required=True, help="Input fragment file name")
    parser.add_argument("-p", "--prefix", dest="prefix", required=True, help="a prefix for output file names")
    parser.add_argument("-f", "--filter", dest="filter", default=None, help="a gene list file that will be used to subset the input data")
    parser.add_argument("-G", "--columnG", dest="columnG", default='gene_id', help="specify the column name for unite feature")
    parser.add_argument("-T", "--columnT", dest="columnT", default='transcript_id', help="specify the column name for feature of which numbers will be counted")
    parser.add_argument("-F", "--columnF", dest="columnF", default='fragment_id', help="specify the column name for feature of which numbers will be counted")
    parser.add_argument("-E", "--columnE", dest="columnE", default='fusion_id', help="specify the column name for feature of which numbers will be counted")
    parser.add_argument("-R", "--columnR", dest="columnR", default='annotation_frequency', help="specify the column name for feature type assigned by EA tool to calculate variable feature ratio, the EA output includes 3 types: Constitutive, Common and Unique. The variable feature ratio is calculated as (Common + Unique) / (Common + Unique + Constitutive)")
    parser.add_argument("-M", "--columnM", dest="columnM", default='flag_multigene', help="specify the column name for multigene flags in feature files, the rows with flags 1 will be removed before the plotting")
    parser.add_argument("-o", "--output", dest="output", required=True, help="Output path")
    parser.add_argument("-s", "--separFeature", dest="separFeature", action='store_false', help="output a separate histogram figure for each feature")
    parser.add_argument("-l", "--lengthhist", dest="lengthhist", action='store_false', help="output a distribution histogram figure for length of each feature")
    parser.add_argument("-r", "--ratiohist", dest="ratiohist", action='store_false', help="output a distribution histogram figure for variable feature ratios for each feature")
    parser.add_argument("-v", "--varFeatureNumber", dest="varFeatureNumber", action='store_false', help="output a distribution histogram figure for variable feature numbers for each feature")
    args = parser.parse_args()
    return(args)


def remove_duprow(col, data):
    n = data[col].tolist()
    temp = []
    dup_index = []
    for i in range(0, len(n)):
        if n[i] not in temp:
            temp.append(n[i])
        else:
            dup_index.append(i)
    data = data.drop(index = dup_index)
    return data


def calculate_number(unit, feature, multigene, data):
    nr = data.shape[0]
    temp = {}
    for i in range(0, nr):
        uni = data[unit][i].split("|")
        fea = data[feature][i].split("|")
        flag = data[multigene][i]
        if flag == 1:
            continue
        for j in uni:
            if j not in temp:
                temp[j] = fea.copy()
            else:
                for k in fea.copy():
                    if k not in temp[j]:
                        temp[j].append(k)
    for key in temp:
        temp[key] = len(temp[key])
    return pd.Series(temp)


def calculate_num_from_pairfile(geneID, transID, annotDF):
    xcrptGene = annotDF[(~annotDF[transID].str.contains("|",regex=False))&
            (~annotDF[geneID].str.contains("|",regex=False))&
            (annotDF[transID]!="Unannotated")][[geneID, transID]].drop_duplicates()
    return xcrptGene.groupby(geneID).size()


def calculate_length(dat, feature):
    start = feature + '_start'
    stop = feature + '_stop'
    ids = feature + "_id"
    res = dat[stop] - dat[start]
    res.index = dat[ids]
    return res

def generate_bins(data, nb = 10, merge_tail = 99.5, ceil = True, hardcut=False):
    bins = [min(data)]
    if not hardcut:
        cut = np.percentile(data, merge_tail)
        itv = (cut - min(data))/(nb-1)
    else:
        cut = 100
        if max(data) > cut:
            itv = (cut - min(data))/(nb-1)
        else:
            cut = np.percentile(data, merge_tail)
            itv = (cut - min(data))/(nb-1)
    
    for i in range(1, nb):
        bins.append(min(data) + i*itv)
    
    if not hardcut:
        bins = bins + [max(data)+0.1]
    else:
        if max(data) > cut:
            bins = bins + [max(data)+0.1]
        else:
            bins = bins + [100]

    if ceil:
        bins = [math.ceil(x) for x in bins]
    return bins


def set_boxAttri(box, color_list):
    for patch, filers, whis, med, color in zip(box['boxes'], box['fliers'], box['whiskers'], box['medians'], color_list):
        patch.set_color(color)
        filers.set_markeredgecolor(color)
        whis.set_color(color)
        med.set_color('black')
    return


def calculate_varFeatureRatio(varList, allList):
    df1 = allList.to_frame()
    df1.columns = ['n']
    df1.insert(0, 'name', allList.index)
    df2 = varList.to_frame()
    df2.columns = ['n']
    df2.insert(0, 'name', varList.index)
    df_merge = df1.merge(df2, how='left', on = 'name', suffixes=('_totalFeature', '_varFeature')).fillna(0)
    df_merge['ratio'] = df_merge['n_varFeature'] / df_merge['n_totalFeature']
    return df_merge

def hist_plot(dat, bins, colors, text=True, omitText = []):
    x_coords = []
    y_coords = []
    first_height = len(dat[dat.values < bins[1]])
    #second_height = len(dat[(dat.values >= bins[1]) & (dat.values < bins[2])])
    #ylim = max(first_height, second_height) * 1.1
    xlim = len(bins) * 1.1
    shrink = (len(bins)-1)/10
    plt.xlim(1, xlim)
    plt.plot(xlim, first_height)
    plt.grid(False)
    plt.xticks([])
    omitIter = 0
    for i in range(0, len(bins)-1):
        x = i+2
        y = len(dat[(dat.values < bins[i+1]) & (dat.values >= bins[i])])
        col = colors[i]
        plt.vlines(x, 0, y, color = col, linewidth=44/shrink)
        if text and (omitIter not in omitText):
            plt.text(x, y, y, fontsize = 10/shrink, horizontalalignment='center',verticalalignment='bottom')
        x_coords.append(x)
        y_coords.append(y)
        omitIter+=1
    plt.ylim(0, max(y_coords)*1.1)
    return x_coords, y_coords



def visualize_feature(ind, out, prefix, title, fname, kbins=10, xlabelrotation=False, firstBin=True, binCeil = True, hardcut = False):
    if firstBin:
        bins = generate_bins(ind[ind>=4] , kbins-2, 99.9, ceil = binCeil, hardcut = hardcut) # kbins (including a merged top 5% bin) + a bin for single transcript genes
        bins = [1, 2, 3] + bins
    else:
        bins = generate_bins(ind, kbins+1, 99.9, ceil = binCeil, hardcut = hardcut)

    fig = plt.figure(figsize=(10,6), facecolor='white',edgecolor='black')
    ax = fig.add_subplot(111)
    x_coords, y_coords = hist_plot(ind, bins, ['steelblue'] * (len(bins)-1))
    tick_pos = [1] + x_coords
    tick_pos = [i+0.5 for i in tick_pos]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([str(round(i, 2)) for i in bins])
    ax.set_title(title, fontsize=14)
    ax.grid(True)
    if xlabelrotation:
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
    outfile = out + '/' + prefix + fname
    plt.savefig(outfile, dpi=600, format='pdf')
    return ax

def counts_projectToColor(lowCol, highCol, counts, firstBinColor = False):
    sumC = sum(counts)
    if sumC == 0:
        return tuple(map(lambda x: x/255 ,list(lowCol)))
    if firstBinColor:
        return tuple(map(lambda x: x/255, list(lowCol)))
    ratio = [i/sumC for i in counts]
    heat_colors = [((lowCol[0] + (highCol[0] - lowCol[0]) * k)/255, (lowCol[1] + (highCol[1] - lowCol[1]) * k)/255, (lowCol[2] + (highCol[2] - lowCol[2]) * k)/255) for k in ratio]
    return heat_colors



def heatmap_ratio_visualize(data, ax, heat_nrows, tick_pos, lowCol=(0, 255, 0), highCol=(255, 0, 0), kb = 10):
    heat_bins = generate_bins([0,1], nb = heat_nrows, merge_tail=100, ceil=False)
    heat_x_coords = heat_bins[0:heat_nrows]
    heat_x_coords = [e+0.05 for e in heat_x_coords]
    heat_starts = tick_pos[0:(kb+1)]
    heat_ends = tick_pos[1:(kb+2)]
    res = []
    for b in range(0, kb+1):
        ind = pd.Series(data[b])
        counts = []
        for i in range(0, len(heat_bins)-1):
            each = ind[(ind.values < heat_bins[i+1]) & (ind.values >= heat_bins[i])]
            counts.append(each.count())
        if b == 0:
            firstBinColor = True
        else:
            firstBinColor = False
        heat_colors = counts_projectToColor(lowCol, highCol, counts, firstBinColor = firstBinColor)
        ax.hlines(heat_x_coords, heat_starts[b], heat_ends[b], color = heat_colors, linewidth = 65/heat_nrows)
        res.append(counts)
    return res


def find_break_lim(y_coords):
    n = len(y_coords)
    sorted_ys = np.flip(np.sort(y_coords))
    sorted_idx = np.flip(np.argsort(y_coords))
    differences = sorted_ys[0:(n-1)] - sorted_ys[1:n]
    max_idx = np.argmax(differences)
    break_y_up = y_coords[sorted_idx[max_idx]]
    break_y_down = y_coords[sorted_idx[max_idx+1]]
    return break_y_up, break_y_down
    
def visualize_combined(list_main, list2, list3, df2_ratio, df3_ratio, out, prefix, kbins=10):
    dat = list_main
    dat2 = list2
    dat3 = list3
    xlabel = 'Number of transcripts per gene'
    ylabel = 'Number of genes'
    ylabel2 = 'Exon regions \n per gene'
    ylabel3 = 'Fragment \n per gene'
    ylabel4 = 'Ratio of variable \n exon region'
    ylabel5 = 'Ratio of \nvariable fragment'
    titile = 'Distribution of number of transcripts per gene'
    #lowColor = (0, 255, 0)
    lowColor = (224, 224, 224)
    highColor = (255, 0, 0)

    font2 = {'family' : 'serif',
             'weight' : 'normal',
             'size'   : 13,
    }
    font3 = {'family' : 'serif',
             'weight' : 'normal',
             'size'   : 10,
    }
    if max(dat) > 18:
        lastbin = max(dat)
    else:
        lastbin = 19
    bins = [1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, lastbin]

    ticklinewidth = 0.5
    colors = sns.color_palette("hls", kbins+1)
    #height = len(dat[(dat.values < bins[1]) & (dat.values >= bins[0])]) * 1.1
    xlim = (kbins + 2) * 1.1
    shrink = (len(bins)-1)/10
    
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(10,10), facecolor='white',edgecolor='black')
    grid = plt.GridSpec(7,1 , wspace=0, hspace=0.15)
    y_coords = []
    x_coords = []
    ###############################
    # plot section 1
    ax1 = fig.add_subplot(grid[0,0])
    ax1.patch.set_facecolor('white')
    ax1.axes.get_yaxis().set_visible(True)
    ax1.spines['bottom'].set_visible(False)
    ax1.set_title(titile, fontsize=14)

    x_coords, y_coords = hist_plot(dat, bins, colors, text=False)

    break_y_up, break_y_down = find_break_lim(y_coords)
    
    ax1.set_ylim(break_y_up * 0.8, max(y_coords) * 1.2)
    for i in range(0, len(bins)-1):
        if y_coords[i] >= break_y_up:
            ax1.text(x_coords[i], y_coords[i], y_coords[i], fontsize = 10/shrink, horizontalalignment='center',verticalalignment='bottom')

    tick_pos = [1] + x_coords
    tick_pos = [i+0.5 for i in tick_pos]

    ##################################
    # plot section 2
    ax2 = fig.add_subplot(grid[1:3,0])
    ax2.patch.set_facecolor('white')
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax2.set_ylabel(ylabel,font2)
    #ax2.plot(x_coords, y_coords, color = 'steelblue', lw = 1) # fit a line
    
    hist_plot(dat, bins, colors, text=False)
    ax2.set_ylim(0, break_y_down * 1.3)
    for i in range(0, len(bins)-1):
        if y_coords[i] <= break_y_down:
            ax2.text(x_coords[i], y_coords[i], y_coords[i], fontsize = 10/shrink, horizontalalignment='center',verticalalignment='bottom')


    ################################
    # extract genes based on bins
    feature2s = []
    feature3s = []
    ratio2s = []
    ratio3s = []
    for i in range(0, len(bins)-1):
        genes = dat[(dat.values < bins[i+1]) & (dat.values >= bins[i])].index.tolist()
        feature2 = dat2[dat2.index.isin(genes)].tolist()
        feature2 = [x for x in feature2 if str(x) != 'nan']
        feature3 = dat3[dat3.index.isin(genes)].tolist()
        feature3 = [x for x in feature3 if str(x) != 'nan']
        feature2s.append(feature2)
        feature3s.append(feature3)
        ratio2 = df2_ratio[df2_ratio['name'].isin(genes)]['ratio'].tolist()
        ratio3 = df3_ratio[df3_ratio['name'].isin(genes)]['ratio'].tolist()
        ratio2s.append(ratio2)
        ratio3s.append(ratio3)

    ################################
    # plot section 3
    ax3 = fig.add_subplot(grid[3,0])
    ax3.patch.set_facecolor('white')
    ax3.spines['top'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.grid(False)
    ax3.set_ylim(0, max(dat2))
    ax3.set_ylabel(ylabel2, font3, labelpad = 11)
    box = ax3.boxplot(feature2s, positions = x_coords, notch = False, patch_artist=True)
    set_boxAttri(box, colors)

    ax3.set_xlim(1, xlim)
    ax3.set_xticklabels([])
    ax3.vlines(tick_pos, 0, max(dat2), color = 'grey', linewidth=ticklinewidth)

    y_arrow = max(dat2)
    for i in range(0, len(bins)-1): 
        ax3.arrow(x_coords[i], y_arrow, 0, -y_arrow/20, overhang=y_arrow/20, head_width=0.2, head_length=1, width = 1, shape="full",fc=colors[i], ec=colors[i],alpha=0.9)


    ##################################
    # plot section 4
    ax4 = fig.add_subplot(grid[4,0])
    ax4.patch.set_facecolor('white')
    ax4.spines['top'].set_visible(False)
    ax4.spines['bottom'].set_visible(False)
    ax4.grid(False)
    ax4.set_ylim(0, max(dat3))
    ax4.set_ylabel(ylabel3, font3, labelpad = 11)
    box2 = ax4.boxplot(feature3s, positions = x_coords, notch = False, patch_artist=True)
    set_boxAttri(box2, colors)

    y_arrow = max(dat3)
    #ax4.set_xticks(tick_pos)
    ax4.set_xticklabels([])
    ax4.set_xlim(1, xlim)
    #ax4.set_xlabel(xlabel, font2, labelpad = 10)
    ax4.vlines(tick_pos, 0, max(dat3), color = 'grey', linewidth=ticklinewidth)
    for i in range(0, len(bins)-1): 
        ax4.arrow(x_coords[i], y_arrow, 0, -y_arrow/20, overhang=y_arrow/20, head_width=0.2, head_length=1, width = 1, shape="full",fc=colors[i], ec=colors[i], alpha=0.9)


    
    ###########################################
    # set heatmap resolution
    heat_nrows = 20
    heat_ylim = 1.1
    
    ####### plot section 5
    ax5 = fig.add_subplot(grid[5,0])
    ax5.patch.set_facecolor('white')
    ax5.spines['top'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax5.grid(False)
    ax5.set_ylim(0, heat_ylim)
    ax5.set_ylabel(ylabel4, font3, labelpad = 11)
    
    hcounts1 = heatmap_ratio_visualize(ratio2s, ax5, heat_nrows, tick_pos, kb = kbins, lowCol = lowColor, highCol = highColor)
    
    ax5.set_xticklabels([])
    ax5.set_xlim(1, xlim)
    ax5.vlines(tick_pos, 0, heat_ylim-0.03, color = '#4682b4', linewidth=ticklinewidth)
    ax5.hlines(0, tick_pos[0], tick_pos[-1], color = '#D3D3D3', linewidth=ticklinewidth*3)

    
     ##################################
    # plot section 6
    ax6 = fig.add_subplot(grid[6,0])
    ax6.patch.set_facecolor('white')
    ax6.spines['top'].set_visible(False)
    ax6.grid(False)
    ax6.set_ylim(0, heat_ylim)
    ax6.set_ylabel(ylabel5, font3, labelpad = 11)
    
    hcounts2 = heatmap_ratio_visualize(ratio3s, ax6, heat_nrows, tick_pos, kb = kbins, lowCol = lowColor, highCol = highColor)

    ax6.set_xlim(1, xlim)
    ax6.set_xlabel(xlabel, font2, labelpad = 10) 
    ax6.set_xticks(tick_pos)
    ax6.set_xticklabels([str(round(i, 2)) for i in bins])
    ax6.vlines(tick_pos, 0, heat_ylim-0.03, color = '#4682b4', linewidth=ticklinewidth)
    
    outfile = out + "/" + prefix + "_combinedfig.pdf"
    outfile2 = out + "/" + prefix + "_combinedfig.png" 
    plt.savefig(outfile, dpi=600, format='pdf')
    plt.savefig(outfile2, dpi=600, format='png')
    return hcounts1, hcounts2, ratio2s, ratio3s


def main():
    args = getoptions()
    dat = pd.read_csv(args.input_g, low_memory=False)
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    if not args.filter == None:
        gene_filter = pd.read_csv(args.filter, header = None, squeeze = True, low_memory=False)
        dat = dat[dat[args.columnG].isin(gene_filter)].reset_index(drop=True)
    num_list_main = calculate_num_from_pairfile(args.columnG, args.columnT, dat)
    del(dat)
    
    dat_e = pd.read_csv(args.input_e, low_memory=False)
    if not args.filter == None:
        dat_e = dat_e[dat_e[args.columnG].isin(gene_filter)].reset_index(drop=True)
    dat_e_var = dat_e[dat_e[args.columnR] != 'Constitutive'].reset_index().drop('index', axis=1)

    num_list_e = calculate_number(args.columnG, args.columnE, args.columnM, dat_e)
    num_list_e_varXcript = calculate_number(args.columnG, args.columnE, args.columnM, dat_e_var)
    len_list_e = calculate_length(dat_e, args.columnE.replace("_id", ""))
    del(dat_e)
    df_ratio_e = calculate_varFeatureRatio(num_list_e_varXcript, num_list_e)
    
    dat_f = pd.read_csv(args.input_f, low_memory=False)
    if not args.filter == None:
        dat_f = dat_f[dat_f[args.columnG].isin(gene_filter)].reset_index(drop=True)
    dat_f_var = dat_f[dat_f[args.columnR] != 'Constitutive'].reset_index().drop('index', axis=1)
    
    num_list_f = calculate_number(args.columnG, args.columnF, args.columnM, dat_f)
    num_list_f_varXcript = calculate_number(args.columnG, args.columnF, args.columnM, dat_f_var)
    len_list_f = calculate_length(dat_f, args.columnF.replace("_id", ""))
    del(dat_f)
    df_ratio_f = calculate_varFeatureRatio(num_list_f_varXcript, num_list_f)
    
    visualize_combined(num_list_main, num_list_e, num_list_f, df_ratio_e, df_ratio_f, out = os.path.normpath(args.output), prefix = args.prefix)
    if args.separFeature:
        visualize_feature(num_list_e, out = os.path.normpath(args.output), prefix = args.prefix, hardcut = True, title = 'Number of exons regions per gene', fname = '_exonRegonsfig.pdf')
        visualize_feature(num_list_f, out = os.path.normpath(args.output), prefix = args.prefix, hardcut = True, title = 'Number of exons fragments per gene', fname = '_fragmentfig.pdf')

    if args.lengthhist:
        visualize_feature(len_list_e, out = os.path.normpath(args.output), firstBin = False, kbins=20, xlabelrotation=True, prefix = args.prefix, title = 'Length distribution of exons regions', fname = '_Length_exonRegion.pdf')
        visualize_feature(len_list_f, out = os.path.normpath(args.output), firstBin = False, kbins=20, xlabelrotation=True, prefix = args.prefix, title = 'Length distribution of exons fragments', fname = '_Length_fragmentfig.pdf')
    
    if args.varFeatureNumber:
        visualize_feature(num_list_f_varXcript, out = os.path.normpath(args.output), kbins=20, prefix = args.prefix, title = 'Number of variable exons fragments per gene', fname = '_variableFragmentfig.pdf')
        visualize_feature(num_list_e_varXcript, out = os.path.normpath(args.output), kbins=20, prefix = args.prefix, title = 'Number of variable exons regions per gene', fname = '_variableExonRegonsfig.pdf')
    
    if args.ratiohist:
        visualize_feature(df_ratio_e['ratio'], out = os.path.normpath(args.output), prefix = args.prefix, firstBin = False, binCeil = False, kbins=20, title = 'Ratio distribution of exons regions', fname = '_RatioVariable_exonRegionfig.pdf')
        visualize_feature(df_ratio_f['ratio'], out = os.path.normpath(args.output), prefix = args.prefix, firstBin = False, binCeil = False, kbins=20, title = 'Ratio distribution of exons fragments', fname = '_RatioVariable_fragmentfig.pdf')
    return


if __name__ == "__main__":
    main()






