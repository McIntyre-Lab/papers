#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def getOptions():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Evaluate expression from flag on off file")

    # Input data
    parser.add_argument("-f", "--flag", dest="inFlag", required=True, help="Input full flag on/off TSV file of novel loci")

    # Output data
    parser.add_argument("-p", "--prefix", dest="outPrefix", required=True, help="Output file prefix")

    args = parser.parse_args()
    return args

def main():
    # Get input flag file
    flagDF = pd.read_csv(args.inFlag,sep="\t")
    
    # Flag genotypes when a least one treatment is detected
    for genotype in ['B73','C123','Hp301','Mo17','NC338']:
        flagDF['flag_'+genotype] = np.where(flagDF['flag_'+genotype+'_Amb']+flagDF['flag_'+genotype+'_Ele']>0,1,0)
    flagDF['sum_genotype_flag'] = flagDF[[c for c in flagDF.columns if ("flag_" in c)&("Amb" not in c)&("Ele" not in c)]].sum(axis=1)
    
    # Get number of genes detected in each group of samples
    countDF = flagDF.groupby([c for c in flagDF.columns if ("flag_" in c)&("Amb" not in c)&("Ele" not in c)])['gene_id'].count().reset_index().rename(columns={'gene_id':'num_gene'}).sort_values('num_gene',ascending=False)
    print("Number of genes detected in at least one treatment of each genotype:\n{}".format(countDF.to_string(index=False)))

    # Plot log10 transformed mean expression values
    logDF = flagDF.set_index('gene_id')[[c for c in flagDF.columns if "mean" in c]].copy()
    for genotype in ['B73','C123','Hp301','Mo17','NC338']:
        logDF['log_mean_'+genotype+"_Amb"] = np.log2(logDF['mean_'+genotype+"_Amb"])
        logDF['log_mean_'+genotype+"_Amb"] = np.where(logDF['log_mean_'+genotype+"_Amb"]==(-np.inf),np.nan,logDF['log_mean_'+genotype+"_Amb"])
        logDF['log_mean_'+genotype+"_Ele"] = np.log2(logDF['mean_'+genotype+"_Ele"])
        logDF['log_mean_'+genotype+"_Ele"] = np.where(logDF['log_mean_'+genotype+"_Ele"]==(-np.inf),np.nan,logDF['log_mean_'+genotype+"_Ele"])
    cols = [c for c in logDF.columns if "log" in c]
    xlabels = list(map(lambda x: "_".join(x.split("_")[2:]),cols))
    logPlot = logDF[cols].plot.box(figsize=(12,12),rot=45)
    logPlot.set_xticklabels(xlabels)
    logPlot.set_ylabel("Log2(mean # reads)")
    plt.show(logPlot)
    plt.tight_layout()
    plt.savefig("{}_log2meanExpression_by_geno_trt.png".format(args.outPrefix),dpi=600,format="png")
    plt.savefig("{}_log2meanExpression_by_geno_trt.pdf".format(args.outPrefix),dpi=600,format="pdf")

    for genotype in ['B73','C123','Hp301','Mo17','NC338']:
        logDF['log_mean_'+genotype] = np.log2(logDF['mean_'+genotype+"_Amb"]+logDF['mean_'+genotype+"_Ele"])
        logDF['log_mean_'+genotype] = np.where(logDF['log_mean_'+genotype]==(-np.inf),np.nan,logDF['log_mean_'+genotype])
    cols = [c for c in logDF.columns if ("log" in c)&("Amb" not in c)&("Ele" not in c)] 
    xlabels = list(map(lambda x: "_".join(x.split("_")[2:]),cols))
    logPlotGeno = logDF[cols].plot.box(figsize=(12,12),rot=45)
    logPlotGeno.set_xticklabels(xlabels)
    logPlotGeno.set_ylabel("Log2(mean # reads)")
    plt.show(logPlotGeno)
    plt.tight_layout()
    plt.savefig("{}_log2meanExpression_by_geno_sum.png".format(args.outPrefix),dpi=600,format="png")
    plt.savefig("{}_log2meanExpression_by_geno_sum.pdf".format(args.outPrefix),dpi=600,format="pdf")

#    meanDF = flagDF.groupby([c for c in flagDF.columns if ("flag_" in c)&("Amb" not in c)&("Ele" not in c)]).agg({
#            'gene_id':'count',
#            'mean_B73_Amb':'describe',
#            'mean_B73_Ele':'describe',
#            'mean_C123_Amb':'describe',
#            'mean_C123_Ele':'describe',
#            'mean_Hp301_Amb':'describe',
#            'mean_Hp301_Ele':'describe',
#            'mean_Mo17_Amb':'describe',
#            'mean_Mo17_Ele':'describe',
#            'mean_NC338_Amb':'describe',
#            'mean_NC338_Ele':'describe'}).reset_index()
    
#    noNC = flagDF[(flagDF['flag_NC338']==0)&(flagDF['sum_genotype_flag']==4)].copy()
#    noNC = noNC.set_index('gene_id')
#    noNC[[c for c in noNC.columns if "mean" in c]].plot.box(figsize=(12,12))
#    onlyNC = flagDF[(flagDF['flag_NC338']==1)&(flagDF['sum_genotype_flag']==1)].copy()
#    onlyNC = onlyNC.set_index('gene_id')
#    onlyNC[[c for c in onlyNC.columns if "mean" in c]].plot.box(figsize=(12,12))
    
#    noB73 = flagDF[(flagDF['flag_B73']==0)&(flagDF['sum_genotype_flag']==4)].copy()
#    noB73 = noB73.set_index('gene_id')
#    noB73[[c for c in noB73.columns if "mean" in c]].plot.box(figsize=(12,12))
    onlyB73 = flagDF[(flagDF['flag_B73']==1)&(flagDF['sum_genotype_flag']==1)].copy()
    onlyB73.to_csv(args.outPrefix+"_B73_only_novel_loci.csv",index=False)
    onlyB73 = onlyB73.set_index('gene_id')
    onlyB73[[c for c in onlyB73.columns if "mean" in c]].plot.box(figsize=(12,12))


if __name__ == '__main__':
    # Parse command line arguments
    global args
    args = getOptions()
    main()
