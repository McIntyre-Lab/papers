#!/bin/sh
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --job-name=fqSplitDups
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/staph_relapse/scripts/SLURM_LOGS/fqSplitDups_%A.%a.out
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=36G
#SBATCH --array=2-79



### run fastqSplitDups on staph relapse data
### S15 sampleNum = Failed Library!!


module load python/2

#Set directories
PROJ=/ufrc/mcintyre/share/staph_relapse
ORIG=${PROJ}/original_data
UNZIP=${PROJ}/unzipped
    mkdir -p ${UNZIP}
FQSPLIT=${PROJ}/fastq_split_dups_PE
    mkdir -p ${FQSPLIT}
LOGS=${FQSPLIT}/logs
    mkdir -p ${LOGS}

## Design file
DESIGN_FILE=${PROJ}/design_files/sampleNum_sampleID_omit_S15_noHeader.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLE_NUM=${ARRAY[0]}
SAMPLEID=${ARRAY[1]}
   ## = 6124-S63_S88_L002_R1_001.fastq.gz ==> sampleNum is S63


## fqSplitDups PE 
echo -e "STARTING FQSPLITDUPS ON PE SAMPLE:  $SAMPLEID"

python ${PROJ}/scripts/fastqSplitDups_2MAI.py \
    -r1 ${UNZIP}/*-${SAMPLE_NUM}_*_R1_001.fastq \
    -r2 ${UNZIP}/*-${SAMPLE_NUM}_*_R2_001.fastq \
    --outdir ${FQSPLIT} \
    -o $LOGS/${SAMPLEID}_fqSplitDup_summary.csv \
    -t $FQSPLIT/${SAMPLEID}_fqSplitDup_table.tsv \
    -g $LOGS/${SAMPLEID}_fqSplitDup.log
