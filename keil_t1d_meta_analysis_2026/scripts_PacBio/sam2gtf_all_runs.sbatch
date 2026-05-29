#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH--job-name=sam2gtf_3rd.%A_%a
#SBATCH--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/sam2gtf_all_runs.%A_%a.out
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --array=2-61
#SBATCH --account=concannon
#SBATCH --qos=concannon-b



### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio


module load samtools/1.10
module load bedtools/2.29.2

## design has header!!
DESIGN_FILE=$PROJ/design_files/dsgn_sampleID_2_bc_2_run.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RUN=${ARRAY[0]} 
RUNID=${ARRAY[1]}
NAME=${ARRAY[2]}
SAMPLEID=${ARRAY[3]}
BC=${ARRAY[4]}

##Set directories

IND=${PROJ}/minimap2_all_runs_hq/${RUN}

OUTD=${PROJ}/sam2gtf_all_runs/
    mkdir -p $OUTD


# convert sam to bam
samtools view -b $IND/${SAMPLEID}_${RUN}_minimap2.sam > $IND/${SAMPLEID}_${RUN}_minimap2.bam

# report each part of split bam alignment as distinct BED interval, create gtf header
bedtools bamtobed \
    -split \
    -i $IND/${SAMPLEID}_${RUN}_minimap2.bam | awk -F "\t" -v sampleID=${SAMPLEID} '{print $1"\t"sampleID"\t""exon""\t"$2+1"\t"$3"\t"".""\t"$6"\t"".""\t""gene_id " "\"" $4 "\"" "; " "transcript_id " "\"" $4 "\"" ";"}' \
    > $OUTD/${SAMPLEID}_${RUN}_minimap2.gtf

rm $IND/${SAMPLEID}_${RUN}_minimap2.sam







