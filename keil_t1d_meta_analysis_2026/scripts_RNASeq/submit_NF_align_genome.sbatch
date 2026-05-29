#!/bin/bash
#SBATCH --job-name=NF_T1DGC_aln
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/blue/concannon/share/T1DGC_R01_controls_RNAseq/scripts/SLURM_LOGS/T1DGC_aln_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=24:00:00
#SBATCH --account=concannon
date;hostname;pwd

## run nextflow script

## directory where nextflow scripts are
NF=/blue/concannon/share/T1DGC_R01_controls_RNAseq/scripts
REPORT=/blue/concannon/share/T1DGC_R01_controls_RNAseq/scripts/nextflow_reports
    mkdir -p ${REPORT}

## create java temp dir
mkdir -p /blue/concannon/share/T1DGC_R01_controls_RNAseq/ROZ_NF

module load nextflow/20.07.1

export _JAVA_OPTIONS="-Djava.io.tmpdir=blue/concannon/share/T1DGC_R01_controls_RNAseq/ROZ_NF"


nextflow run ${NF}/NF_align_genome.nf \
    -c ${NF}/nextflow_align_genome.config \
    -with-report ${REPORT}/report_NF_align_genome.html
