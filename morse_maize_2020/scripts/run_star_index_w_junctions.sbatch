#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=star
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/star-%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32gb
#SBATCH --cpus-per-task=8
#SBATCH --account=mcintyre


### star to generate junctions
	## readLength is 150bp
	## length on either side of junction is:  (RL/2 + 13)
        ## juncLength is:  (150/2 + 13) *2 = 176
    module load star



    ### Set Directories
    PROJ=/ufrc/mcintyre/share/maize_ainsworth
    NUMPROC=8
    REFDIR=$PROJ/pacbio_analysis/references/star_index
        mkdir -p $REFDIR

    GTF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/gtf/zea_mays/Zea_mays.B73_RefGen_v4.41.gtf
    FASTA=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa
    JUNCLENGTH=176
    
    ### trying readLength-1 for sjdbOverhang
   
#    STAR --runThreadN ${NUMPROC} \
#        --runMode genomeGenerate \
#        --genomeDir ${REFDIR} \
#        --genomeFastaFiles ${FASTA} \
#        --sjdbGTFfile ${GTF} \
#        --sjdbGTFtagExonParentTranscript Parent \
#        --sjdbOverhang 149

    STAR --runThreadN ${NUMPROC} \
        --runMode genomeGenerate \
        --genomeDir ${REFDIR} \
        --genomeFastaFiles ${FASTA} \
        --sjdbGTFfile ${GTF}

