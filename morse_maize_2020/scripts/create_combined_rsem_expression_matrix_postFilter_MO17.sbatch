#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=matrixC
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/matrixC_%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=4
#SBATCH --account=mcintyre


### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
SCRIPTS=$PROJ/scripts

OUTDIR=$PROJ/rsem_output_postFilter_MO17/isoform_results
  
ROZ=$PROJ/rozEM
    mkdir -p $ROZ

#### merging rsem bowtie1 isoform.results files -- expected_counts
	## pull out only transcriptID and expected_counts columns from each file
	## rename expected_counts column to the sample name
	## paste all files together and keep transcriptID and expected counts for each
    	## run script that combines SE and PE expression counts 


    for i in SE PE
    do
        for samp in $OUTDIR/*_${i}.isoforms.results
        do
          name=$(basename ${samp})
          echo ${name}
          cut -f1,5 ${samp} > $ROZ/${name}".expected_count"
          s=$(echo "${name}" | sed 's/_'"${i}"'.isoforms.results//')
          echo ${s}
          sed s/'expected_count'/${s}/g $ROZ/${name}".expected_count" > $ROZ/${s}_${i}.exp
        done

	## paste all output together and keep transcriptID and expected count for each sample
        paste $ROZ/*_${i}.exp > $ROZ/${i}_all.exp

        ## cols variable is a comma separated list of all columns to keep in $ROZ/${i}_all.exp
        COLS=$(awk '{ cols = "1"; \
            for (x=1; x<=NF; x++) \
            {if (x % 2 ==0) \
                { cols = cols","x; \
                } } } END {print cols}' $ROZ/${i}_all.exp)
        echo $COLS

        cut -f $COLS $ROZ/${i}_all.exp > $ROZ/${i}_expression_matrix.post-filter.expected_count.txt
    done

    module purge 
    module load python3
  
    ## run script that combines (sums) SE and PE expected counts  
    python3 $SCRIPTS/merge_rsem_expression.py \
        -p $ROZ/PE_expression_matrix.post-filter.expected_count.txt \
        -s $ROZ/SE_expression_matrix.post-filter.expected_count.txt \
        -o $OUTDIR/combined_expression_matrix.post-filter.expected_count.txt

rm -r $ROZ
