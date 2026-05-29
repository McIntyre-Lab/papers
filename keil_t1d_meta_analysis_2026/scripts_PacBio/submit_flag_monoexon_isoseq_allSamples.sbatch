#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=flag_mono_isoseq2
#SBATCH --output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/flag_mono_isoseq2.%j.out
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1

module load python/3.6

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
SCRIPTS=$PROJ/scripts
OUTD=$PROJ/flag_monoexon_isoseq_allSamples_02ndk
   mkdir -p ${OUTD}
IND=${PROJ}/corrected_gtf_isoseq_allSamples_02ndk

REF=/blue/mcintyre/share/transcript_distance/human_analysis/TranD_pairwise_ensembl_refseq_union/consolidated_transcriptome.gtf

 ## Flag monoexon reads and references
    python ${SCRIPTS}/flag_monoexon_03avn.py \
        -r ${REF} \
        -i ${IND}/isoseq_allSamples_corrected_geneID.gtf \
        -o ${OUTD}/flag_monoexon_isoseq_allSamples.csv
