#!/bin/bash

#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=blast
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/blast_%A.out
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b

### blast non-b73 and non-mo17 samples that do NOT map to either B73 nor to mo17 to ncbi nr database


module load ncbi_blast python

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTPUT=$PROJ/evidence_genotype_specific_seqs

ROZ=$PROJ/rozBlast
    mkdir -p $ROZ


### query = unmapped fasta 
### subject = nr database

for FA in $OUTPUT/seqs_no_aln_*.fa
do
        SAMPLE=$(basename $FA)
        echo ${FA}
        echo ${SAMPLE}

        QUERY=${FA}
        
        ## create header for blast tsv output
        echo "qseqid	stitle	qstart	qend	qlen	sstart	send	slen	length	bitscore	evalue	mismatch	gaps" > $ROZ/header.tsv

        ## blast and output csv file
        blastn \
            -db nt \
            -query $QUERY \
            -outfmt "6 qseqid stitle qstart qend qlen sstart send slen length bitscore evalue mismatch gaps" \
            -out $ROZ/blast_${SAMPLE}_vs_nt.tsv

        ## add header line
            cat $ROZ/header.tsv $ROZ/blast_${SAMPLE}_vs_nt.tsv > $OUTPUT/blast_${SAMPLE}_vs_nt.tsv

        ### blast and output html file
            blastn \
            -db nt \
            -query $QUERY \
            -html \
            -out $OUTPUT/blast_${SAMPLE}_vs_nt.html
done

rm -r $ROZ
