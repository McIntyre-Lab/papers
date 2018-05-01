# Load packages

import argparse
import numpy
import matplotlib
import pandas
import math
import os
import scipy
from scipy import stats, integrate

import seaborn

seaborn.set(color_codes=True)
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Polygon

parser = argparse.ArgumentParser(description="Import CSV and give output names")
parser.add_argument("--input", dest="csvInput", action='store', required=True, help="Input CSV")
parser.add_argument("--output-corr", dest="RepOut", action='store', required=True, help="Output filename")
parser.add_argument("--output-ba", dest="BAOut", action='store', required=True, help="Output filename")
args = parser.parse_args()


#import data

mcLab=os.environ["MCLAB"]


rsemData=numpy.genfromtxt(args.csvInput, delimiter=",", names=True,
                        dtype=[('transcript_id','|S30'),('tpm_nsc1','<f8'),('tpm_nsc2','<f8'),('log_tpm_nsc1','<f8'),      
                               ('log_tpm_nsc2','<f8'),('mean_tpm','<f8'),('min_log_tpm','<f8'),('mean_log_tpm','<f8'),
                               ('min_tpm','<f8'),('tpm_bin','<f8'),('diff_tpm','<f8'),('diff_log_tpm','<f8'),
                               ('diff_over_mean_log_tpm','<f8')])
#set up data
logRep1=[]
logRep2=[]
diffLog=[]
binColor=[]
meanLog=[]

## Trying bins : min<0.5, min<1, min<2, other

for i in range(0,len(rsemData)):
    logRep1.append(rsemData[i][3])
    logRep2.append(rsemData[i][4])
    diffTPM=rsemData[i][3]-rsemData[i][4]
    diffLog.append(diffTPM)
    meanLog.append(rsemData[i][7])
    if rsemData[i][9] == 0 :
        binColor.append("#000000")
    elif rsemData[i][9] == 1 :
        binColor.append("#e82222")
    elif rsemData[i][9] == 2 :
        binColor.append("#F57900")
    elif rsemData[i][9] == 3 :
        binColor.append("#f4bf42")
    elif rsemData[i][9] == 4 :
        binColor.append("#2274e8")
    else :
        binColor.append("#ff00ff")


fig,ax = plt.subplots(figsize=(6,4), dpi=300)
pltCorr2 = ax.scatter(logRep1,logRep2, s=20, color=binColor)

plt.xlim(-0.5,10.5)
plt.ylim(-0.5,10.5)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('log-TPM NPC Replicate 2')
ax.set_xlabel('log-TPM NPC Replicate 1') 
ax.set_title('Concordance of RSEM estimation between NPC replicates')

fig.savefig(args.RepOut)


fig1,ax = plt.subplots(figsize=(6,4), dpi=300)
pltCorr2 = ax.scatter(meanLog,diffLog, s=20, color=binColor)
plt.xlim(-0.5,10.5)
plt.ylim(-6.5,6.5)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('logTPM NPC Rep 1 - logTPM NPC Rep 2')
ax.set_xlabel('mean logTPM (NPC Rep 1, TPM NPC Rep 2)') 
ax.set_title('Bland-Altman plot')


fig1.savefig(args.BAOut)

