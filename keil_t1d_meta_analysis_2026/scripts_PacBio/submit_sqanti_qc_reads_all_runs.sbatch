#!/bin/bash
#SBATCH --mail-user=nkeil@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=sqanti3QC_ENS
#SBATCH --output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/sqanti3QC_reads.%A-%a.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH --cpus-per-task=12
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-61

module load sqanti3/4.2

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio


## design has header!!
DESIGN_FILE=$PROJ/design_files/dsgn_sampleID_2_bc_2_run.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

RUN=${ARRAY[0]} 
RUNID=${ARRAY[1]}
NAME=${ARRAY[2]}
SAMPLEID=${ARRAY[3]}
BC=${ARRAY[4]}


#Set directories and files

UNION_GTF=/blue/mcintyre/share/transcript_distance/human_analysis/TranD_pairwise_ensembl_refseq_union/consolidated_transcriptome.gtf
REF_FA=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa

IND=${PROJ}/sam2gtf_all_runs
OUTD=${PROJ}/sqanti3_qc_reads_all_runs
    mkdir -p ${OUTD}

##  --gtf option indicates GTF isoform input
##  --skipORF option to skip ORF prediction (runs much faster)
##  --force_id_ignore option to allow transcript ID names that are not like PB.X.Y
##  Order of files is (1) isoforms (FA/FQ/GTF), (2) annotation (ref GTF/GFF), (3) genome (ref FA)

#   --dir is path to output directory
#   -o prefix name to use

## order of files input into sqanti3_qc.py is critical!!
# 1. sam2gtf (from read align to genome)
# 2. FSM consolidated ref gtf file 
# 3. ref genome fasta file 


sqanti3_qc.py \
    --skipORF \
    --force_id_ignore \
    --report pdf \
    --dir ${OUTD} \
    -o ${SAMPLEID}_${RUN} \
    ${IND}/${SAMPLEID}_${RUN}_minimap2.gtf \
    ${UNION_GTF} \
    ${REF_FA}


