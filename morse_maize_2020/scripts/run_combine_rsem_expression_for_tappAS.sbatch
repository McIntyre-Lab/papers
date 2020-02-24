#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL
#SBATCH --job-name=combMatrix
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/combMatrix_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre
#SBATCH --dependency=afterok:41507491

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
INPUT=$PROJ/RNAseq_rsem_expression_subset_fsm_ism_nic_nnc
SCRIPTS=/ufrc/mcintyre/share/maize_ainsworth/scripts

ROZ=$INPUT/roz_combine_rsem
    mkdir -p ${ROZ}

#### merging rsem bowtie1 genes.results and isoforms.results files -- expected_counts and TPM
	## pull out only geneID and expected_counts columns from each file
	## rename expected_counts column to the sample name
	## paste all files together and keep transcriptID and expected counts for each
    	## run script that combines SE and PE expression counts 

    for TYPE in isoforms genes
    do
#        for i in SE PE
#        do
#            for samp in $INPUT/*_${i}.${TYPE}.results
#            do
#            ### expected counts
#                name=$(basename ${samp})
#                echo ${name}
#                cut -f1,5 ${samp} > ${ROZ}/${name}".expected_count"
#                s=$(echo "${name}" | sed s/_${i}.${TYPE}.results//)
#                echo ${s}
#                sed s/'expected_count'/${s}/g ${ROZ}/${name}".expected_count" > ${ROZ}/${s}_${i}.exp
#            done

            ## paste all output together and keep transcriptID and expected count for each sample
#            paste ${ROZ}/*_${i}.exp > $ROZ/${i}_all.exp

            ## cols variable is a comma separated list of all columns to keep in $ROZ/${i}_all.exp
#            COLS=$(awk '{ cols = "1"; \
#                for (x=1; x<=NF; x++) \
#                {if (x % 2 ==0) \
#                    { cols = cols","x; \
#                    } } } END {print cols}' $ROZ/${i}_all.exp)
#            echo ${COLS}

#            cut -f ${COLS} ${ROZ}/${i}_all.exp > ${INPUT}/${i}_${TYPE}_expression_matrix.filtered_GTF.expected_count.txt

#            for samp in ${INPUT}/*_${i}.${TYPE}.results
#            do
            ### TPM
#                name=$(basename ${samp})
#                echo ${name}
#                cut -f1,6 ${samp} > $ROZ/${name}".TPM"
#                s=$(echo "${name}" | sed s/_${i}.${TYPE}.results//)
#                echo ${s}
#                sed s/'TPM'/${s}/g $ROZ/${name}".TPM" > ${ROZ}/${s}_${i}.TPM
#            done

            ## paste all output together and keep transcriptID and expected count for each sample
#            paste ${ROZ}/*_${i}.TPM > ${ROZ}/${i}_all.TPM

            ## cols variable is a comma separated list of all columns to keep in $ROZ/${i}_all.TPM
#            COLS=$(awk '{ cols = "1"; \
#                for (x=1; x<=NF; x++) \
#                {if (x % 2 ==0) \
#                    { cols = cols","x; \
#                    } } } END {print cols}' $ROZ/${i}_all.TPM)
#            echo ${COLS}

#            cut -f ${COLS} ${ROZ}/${i}_all.TPM > ${INPUT}/${i}_${TYPE}_expression_matrix.filtered_GTF.TPM.txt 

#            rm $ROZ/*

#        done

        module purge 
        module load python3
  
        ## run script that combines (sums) SE and PE expected counts  
        python3 ${SCRIPTS}/merge_rsem_expression_${TYPE}.py \
            -p ${INPUT}/PE_${TYPE}_expression_matrix.filtered_GTF.expected_count.txt \
            -s ${INPUT}/SE_${TYPE}_expression_matrix.filtered_GTF.expected_count.txt \
            -o ${INPUT}/combined_expression_${TYPE}_matrix.filtered_GTF.expected_count.txt

        ## run script that combines (sums) SE and PE TPM
        python3 ${SCRIPTS}/merge_rsem_expression_${TYPE}.py \
            -p ${INPUT}/PE_${TYPE}_expression_matrix.filtered_GTF.TPM.txt \
            -s ${INPUT}/SE_${TYPE}_expression_matrix.filtered_GTF.TPM.txt \
            -o ${INPUT}/combined_expression_${TYPE}_matrix.filtered_GTF.TPM.txt

#        rm ${INPUT}/PE_${TYPE}_*
#        rm ${INPUT}/SE_${TYPE}_*
    done

    ### create file containing isoform lengths
    cut -f1-4 ${INPUT}/B73_P1_C1_Amb_SE.isoforms.results > ${INPUT}/filtered_PB.isoforms_lengths.txt


#rm -r $ROZ
