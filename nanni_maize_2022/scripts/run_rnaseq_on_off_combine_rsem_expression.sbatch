#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=on_off
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/rnaseq_on_off_PB_rsem_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --dependency=afterok:41508641

module purge
module load python

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
INPUT=$PROJ/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc

#### merging rsem STAR genes.results and isoforms.results files -- expected_counts and TPM
	## pull out only geneID and expected_counts columns from each file
	## rename expected_counts column to the sample name
	## paste all files together and keep transcriptID and expected counts for each
    	## run script that combines SE and PE expression counts 

    for TYPE in genes isoforms; do
        if [[ ${TYPE} == "isoforms" ]]; then
            NAME="transcript"
        else
            NAME="gene"
        fi
        for RESULT in expected_count TPM; do
            ## Set minimum values
            if [[ ${RESULT} == "expected_count" ]]; then
                MIN=5
            else
                MIN=1
            fi

            ## Get on/off flags based on minimum values set above
            echo ${INPUT}/combined_expression_${TYPE}_matrix.filtered_GTF.${RESULT}.txt
            python $PROJ/scripts/flag_on_off_rsem.py \
                -i ${INPUT}/combined_expression_${TYPE}_matrix.filtered_GTF.${RESULT}.txt \
                -m ${MIN} \
                -t ${NAME} \
                -s maize \
                -o ${INPUT}/flag_on_off_rsem_expression_${TYPE}.filtered_GTF.${RESULT}
        done
    done
done
