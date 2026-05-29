#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH	--job-name=sqanti3QC
#SBATCH	--output=/blue/concannon/share/T1D_R01_controls_PacBio/scripts/SLURM_LOGS/sqanti3QC_%A-%a.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=12gb
#SBATCH --cpus-per-task=12
#SBATCH --account=concannon
#SBATCH --qos=concannon-b
#SBATCH --array=2-13

## convert ref fsm consolidated gff to gtf
## run sqanti3 QC

module load gffread/0.11.8 sqanti3/1.6

### Set Directories
PROJ=/blue/concannon/share/T1D_R01_controls_PacBio
SAM2GTF=${PROJ}/sam2gtf
FSM=${PROJ}/FSM_consolidation

REF_FA=/blue/concannon/share/jnewman/references/hg38_ensembl/hg38.nochr.fa

OUTD=${PROJ}/sqanti3_qc_allRuns
    mkdir -p $OUTD

## convert consolidated	FSM GFF	to GTF format
#sort -k1,1 -k3,3 -k4,4n -k5,5nr ${FSM}/hg38_150bp_FSM_consolidation_all_transcript.gff | gffread - -T -o ${FSM}/hg38_150bp_FSM_consolidation_all_transcript.gtf  

REF_GTF=${FSM}/hg38_150bp_FSM_consolidation_all_transcript.gtf

## design has header!!
DESIGN_FILE=$PROJ/design_files/design_sampleID_2_barcode.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

SAMPLEID=${ARRAY[0]}
BC=${ARRAY[2]}


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


for RUN in first second third fourth fifth
do 

PREFIX=${SAMPLEID}_${RUN}

if [[ ${RUN} == 'first' ]]
then
sqanti3_qc.py \
    --gtf \
    --skipORF \
    --force_id_ignore \
    --dir ${OUTD} \
    -o ${PREFIX} \
    ${SAM2GTF}/${SAMPLEID}_minimap2_${BC}.gtf \
    ${REF_GTF} \
    ${REF_FA}
else
sqanti3_qc.py \
    --gtf \
    --skipORF \
    --force_id_ignore \
    --dir ${OUTD} \
    -o ${PREFIX} \
    ${SAM2GTF}/${RUN}_run/${SAMPLEID}_minimap2_*.gtf \
    ${REF_GTF} \
    ${REF_FA}
fi

done


