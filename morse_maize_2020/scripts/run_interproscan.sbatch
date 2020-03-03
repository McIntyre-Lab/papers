#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=interproscan
#SBATCH --output=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB/scripts/SLURM_LOGS/interproscan_%A_%a.out
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --qos=mcintyre
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-20

module load java/1.8.0_31 python perl

## Require the interproscan module to run interproscan.sh
module load iprscan/5.29-68.0

### InterproScan of filtered SQANTI QC .faa output files
###    output files were split into 20 different files for analysis

### Set Directories
PROJ=/ufrc/mcintyre_conesa_transvar/isoAnnot/maize_ainsworth_PB
IN=$PROJ/split_faa
OUT=$PROJ/interproscan
    mkdir -p ${OUT}

## Split array jobs by files in split_faa
DESIGN_FILE=$PROJ/design_files/df_split_faa.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

FAA=${ARRAY[0]}

## Move into output directory to ensure the temp and tmp
## 	directories created in these steps are not randomly placed
cd ${OUT}


## Run InterProScan
## Added TSV output format to match what is in instructions4tappas.cfg
## Removed TMHMM-2.0c and SignalP_EUK-4.1 applications from options
## 	these require separate license to run in interproscan
##	so they will be run separately
interproscan.sh -i ${IN}/${FAA} -d ${OUT} --disable-precalc -f TSV,XML,HTML \
	-goterms -iprlookup -pa -appl Coils-2.2.1,Pfam-31.0,MobiDBLite-1.0 


## Waited to load these modules because they require older versions
## 	of dependencies
module load trinotate/3.0.1 signalp/4.1

## Run TMHMM
## 	in the trinotate module
## 	Outputs a folder named TMHMM_##### (random numbers at the end)
## 		which contains output .plp and .gnuplot files
TM=${OUT}/${FAA}.tmhmm.out
tmhmm --workdir ${OUT} < ${IN}/${FAA} > ${TM}

## Parse TMHMM
TMPARSE=${OUT}/${FAA}.tmhmm.parsed.out
awk '{if($1!="#" && $3=="TMhelix"){print $1"\tTMHMM\tDOMAIN\tTMhelix\tTMhelix\tRegion of a membrane-bound protein predicted to be embedded in the membrane.\t"$4"\t"$5"\tNA\tNA"}}'\
	${TM} > ${TMPARSE}

## Run SignalP
## Use summary format to retreive D values
## 	(aka positions of possible signal peptide)
SIGP=${OUT}/${FAA}.euk10.summary.out
signalp -t euk -f summary ${IN}/${FAA} > ${SIGP}

## Parse SignalP Summary
## Select only those that were scored to have signal peptides
## 	(aka 'YES' under 'signal peptide?' column)
SIGPparse=${OUT}/${FAA}.euk10.parsed.out
awk '{split($0,a,"="); if($1=="D"){score=$0;} if(a[1]=="Name")\
	{print $1"\t"$NF"\t"score}}' ${SIGP} | awk '{if($NF=="YES")\
	{split($1,name,"="); split($2,networks,"="); split($4,pos,"-");\
	print name[2]"\tSIGNAL_EUK\tDOMAIN\t"networks[2]"\t"networks[2]\
	"\tNone\t"pos[1]"\t"pos[2]"\tNA\tNA"}}' > ${SIGPparse}

## Remove temp and tmp directories made during analyses if empty
if [ -d "tmp" ]; then
	if [ -z "$(ls -A tmp)" ]; then
		rm -r tmp
	fi
fi

if [ -d	"temp" ]; then
        if [ -z	"$(ls -A temp)" ]; then
                rm -r temp
        fi
fi
