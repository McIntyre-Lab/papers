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

### blast mo17 samples that do NOT map to B73 but DO to mo17 yan to b73 genome


module load ncbi_blast python

### Set Directories
PROJ=/ufrc/mcintyre/share/maize_ainsworth
OUTPUT=$PROJ/evidence_genotype_specific_seqs
REF=/ufrc/mcintyre/share/references/maize_b73/b73_genome_ensemblgenomes/fasta/zea_mays/dna

ROZ=$PROJ/rozBlast
    mkdir -p $ROZ


### query = unmapped fasta 
### subject = B73 genome

SUBJECT=$REF/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.BLAST

for FA in $OUTPUT/seqs_yes_aln_*.fa
do
        SAMPLE=$(basename $FA)
        echo ${FA}
        echo ${SAMPLE}

        QUERY=${FA}
        
        ## create header for blast tsv output
        echo "qseqid	stitle	qstart	qend	qlen	sstart	send	slen	length	bitscore	evalue	mismatch	gaps" > $ROZ/header.tsv

        ## blast and output csv file
        blastn \
            -db $SUBJECT \
            -query $QUERY \
            -outfmt "6 qseqid stitle qstart qend qlen sstart send slen length bitscore evalue mismatch gaps" \
            -out $ROZ/blast_${SAMPLE}_vs_b73_genome.tsv

        ## add header line
            cat $ROZ/header.tsv $ROZ/blast_${SAMPLE}_vs_b73_genome.tsv > $OUTPUT/blast_${SAMPLE}_vs_b73_genome.tsv

        ### blast and output html file
            blastn \
            -db $SUBJECT \
            -query $QUERY \
            -html \
            -out $OUTPUT/blast_${SAMPLE}_vs_b73_genome.html
done

rm -r $ROZ
