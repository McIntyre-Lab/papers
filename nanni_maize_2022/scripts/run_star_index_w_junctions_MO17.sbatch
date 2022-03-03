#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=starJunc
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/starJunc-%A.out
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem=56gb
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
    REFDIR=$PROJ/pacbio_analysis/references/star_index_mo17
        mkdir -p $REFDIR

    GTF=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN/Zm-Mo17-REFERENCE-YAN-1.0.mw2.gtf
    FASTA=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN/Zm-Mo17-REFERENCE-YAN-1.0.fasta
    
    STAR --runThreadN ${NUMPROC} \
        --runMode genomeGenerate \
        --genomeDir ${REFDIR} \
        --genomeFastaFiles ${FASTA} \
        --sjdbGTFfile ${GTF}

