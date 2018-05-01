#!/bin/sh
#SBATCH --mail-user=jrbnewman@ufl.edu
#SBATCH --partition=hpg1-compute
#SBATCH --account=mcintyre
#SBATCH --job-name=rsem_refs
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /ufrc/mcintyre/share/conesa_isoform_check/scripts/SLURM_LOGS/out.%j.%A.%a.out
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12gb
#SBATCH --array=5-8
#

## 7

module load R
module load bowtie/0.12.9
module load rsem/1.2.28

PROJ=/ufrc/mcintyre/share/conesa_isoform_check
REF=$PROJ/references/rsem



    DESIGN_FILE=$PROJ/design_files/rsem_references2.csv
    DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    REFSET=${ARRAY[0]}
    #GENE2XS=${REFSET}_xs2gene.txt
    GENE2XS=all_xs2gene.txt
    FASTA=${REFSET}.fa
    RSEMREF=${REFSET}_RSEM

rsem-prepare-reference --bowtie -p 2 \
                       --transcript-to-gene-map $REF/${GENE2XS} \
                       $REF/${FASTA} \
                       $REF/${RSEMREF}

