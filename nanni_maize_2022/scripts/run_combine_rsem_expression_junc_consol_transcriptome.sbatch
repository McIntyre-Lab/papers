#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=rsem_calc
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio_transcriptome_eval/SLURM_LOGS/combine_expression_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre

## Combine individual expression files into expression matrices
##     for expected_count and TPM (isoform and gene level)

module purge
module load python/3.7

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
INDIR=$PROJ/RNAseq_rsem_fsm_ism_nic_nnc_consol_junc
ROZ=${INDIR}/roz_combine_rsem
    mkdir -p ${ROZ}

#### merging rsem star genes.results and isoforms.results files -- expected_counts and TPM
        ## pull out only geneID and expected_counts columns from each file
        ## rename expected_counts column to the sample name
        ## paste all files together and keep transcriptID and expected counts for each
        ## run script that combines SE and PE expression counts


for TYPE in isoforms genes; do
    for i in SE PE; do
        ## expected_count
        for SAMP in ${INDIR}/*_${i}.${TYPE}.results; do
            NAME=$(basename ${SAMP})
            echo ${NAME}
            cut -f 1,5 ${SAMP} > ${ROZ}/${NAME}.expected_count
            s=$(echo ${NAME} | sed s/_${i}.${TYPE}.results//)
            echo ${s}
            sed s/'expected_count'/${s}/g ${ROZ}/${NAME}.expected_count > ${ROZ}/${s}_${i}.exp
        done
        ## paste all output together and keep transcriptID and expected count for each sample
        paste ${ROZ}/*_${i}.exp > ${ROZ}/${i}_all.exp

        ## cols variable is a comma separated list of all columns to keep in $ROZ/${i}_all.exp
        COLS=$(awk '{ cols = "1"; \
            for (x=1; x<=NF; x++) \
            {if (x % 2 ==0) \
                { cols = cols","x; \
                } } } END {print cols}' $ROZ/${i}_all.exp)
        echo ${COLS}

        cut -f ${COLS} ${ROZ}/${i}_all.exp > ${INDIR}/${i}_${TYPE}_expression_matrix.expected_count.txt

        ## TMP
        for SAMP in ${INDIR}/*_${i}.${TYPE}.results; do
            NAME=$(basename ${SAMP})
            echo ${NAME}
            cut -f 1,6 ${SAMP} > ${ROZ}/${NAME}.TPM
            s=$(echo ${NAME} | sed s/_${i}.${TYPE}.results//)
            echo ${s}
            sed s/'TPM'/${s}/g ${ROZ}/${NAME}.TPM > ${ROZ}/${s}_${i}.TPM
        done
        ## paste all output together and keep transcriptID and TPM for each sample
        paste ${ROZ}/*_${i}.TPM > ${ROZ}/${i}_all.TPM

        ## cols variable is a comma separated list of all columns to keep in $ROZ/${i}_all.TPM
        COLS=$(awk '{ cols = "1"; \
            for (x=1; x<=NF; x++) \
            {if (x % 2 ==0) \
                { cols = cols","x; \
                } } } END {print cols}' $ROZ/${i}_all.TPM)
        echo ${COLS}

        cut -f ${COLS} ${ROZ}/${i}_all.TPM > ${INDIR}/${i}_${TYPE}_expression_matrix.TPM.txt

    done

    ## run script that combines (sums) SE and PE expected counts
    python $PROJ/scripts/merge_rsem_expression_${TYPE}.py \
        -p ${INDIR}/PE_${TYPE}_expression_matrix.expected_count.txt \
        -s ${INDIR}/SE_${TYPE}_expression_matrix.expected_count.txt \
        -o ${INDIR}/combined_${TYPE}_expression_matrix.expected_count.txt

    ## run script that combines (sums) SE and PE TPM
    python $PROJ/scripts/merge_rsem_expression_${TYPE}.py \
        -p ${INDIR}/PE_${TYPE}_expression_matrix.TPM.txt \
        -s ${INDIR}/SE_${TYPE}_expression_matrix.TPM.txt \
        -o ${INDIR}/combined_${TYPE}_expression_matrix.TPM.txt

    rm ${INDIR}/PE_${TYPE}_*
    rm ${INDIR}/SE_${TYPE}_*
done

### create file containing isoform lengths
cut -f 1-4 ${INDIR}/B73_P1_C1_Amb_SE.isoforms.results > ${INDIR}/fsm_ism_nic_nnc_consol_junc.isoforms_lengths.txt

rm -r ${ROZ}
