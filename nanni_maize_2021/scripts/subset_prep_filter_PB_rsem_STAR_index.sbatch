#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=rsem_index
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/pacbio_rsem_STAR_index_%A.out
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b

module purge
module load star/2.7.0b rsem/1.2.31

## Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUT=$PROJ/pacbio_analysis/subset_FSM_ISM_NIC_NNC/subset_references
    mkdir -p ${OUT}

FA=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa

### testing if sqanti gtf will index without subsettin
###GTF=/ufrc/mcintyre/share/maize_ainsworth/sqanti_post_filter/sqanti_filtered_corrected.gtf
GTF=$PROJ/pacbio_analysis/subset_FSM_ISM_NIC_NNC/pbID_fsm_ism_nic_nnc.gtf


## Make rsem STAR index
echo -e "Prepare rsem STAR index... $(date)\n"

rsem-prepare-reference --star -p 8 \
    --gtf ${GTF} ${FA} $OUT/subset_fsm_ism_nic_nnc_rsem_index_STAR
   
echo "prep STAR ref done."
