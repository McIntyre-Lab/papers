#!/bin/bash
#SBATCH --job-name=nf
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/rnaseq/SLURM_LOGS/nf_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --time=12:00:00
#SBATCH --account=mcintyre
date;hostname;pwd

## run nextflow script
## directory where nextflow scripts are
NF=/blue/mcintyre/share/maize_ainsworth/scripts/rnaseq
REPORT=/blue/mcintyre/share/maize_ainsworth/mo17_align_nextflow_reports
    mkdir -p ${REPORT}

## create java temp dir
mkdir -p /blue/mcintyre/share/maize_ainsworth/ROZ_NF_mo17_cau

module load nextflow/20.07.1

export _JAVA_OPTIONS="-Djava.io.tmpdir=/blue/mcintyre/share/maize_ainsworth/ROZ_NF_mo17_cau"

nextflow run ${NF}/NF_rnaseq_align_genome_coverage_maize.nf \
    -c ${NF}/nextflow_run_rnaseq_align_genome_coverage_mo17_cau.config \
    -with-report ${REPORT}/align_genome_coverage_mo17_cau.html

