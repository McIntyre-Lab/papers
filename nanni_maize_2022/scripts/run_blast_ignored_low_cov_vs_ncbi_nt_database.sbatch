#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=blast
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/blast_cov_%A_%a.out
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --array=1-11

### blast PB transcripts low coverage to b73, mo17 yand or mo17 cau references to ncbi nr database

module load ncbi_blast python

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTPUT=$PROJ/compare_b73_2_mo17/blast_low_coverage_PBreads
    mkdir -p $OUTPUT

ROZ=$PROJ/rozBC
    mkdir -p $ROZ

for REF in b73 mo17_yan mo17_cau
#for REF in mo17_yan mo17_cau
do

    ## Design file
    DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
    DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
    IFS=',' read -ra ARRAY <<< "$DESIGN"

    GENO=${ARRAY[1]}
    TRT=${ARRAY[2]}    
    ID=${ARRAY[3]}

    SAMPLE=${ID}_${GENO}_${TRT}

    ### query = low coverage fasta 
    ### subject = nr database

    FA=$PROJ/compare_b73_2_mo17/uncollapsed_${REF}_IsoSeq_reads/${SAMPLE}_unCollapsed_low_cov_${REF}_ref.fa

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
        -out $ROZ/blast_low_cov_${SAMPLE}_2_${REF}_vs_nt.tsv

    ## add header line
    cat $ROZ/header.tsv $ROZ/blast_low_cov_${SAMPLE}_2_${REF}_vs_nt.tsv > $OUTPUT/blast_low_cov_${SAMPLE}_2_${REF}_vs_nt.tsv

    ### blast and output html file
    blastn \
        -db nt \
        -query $QUERY \
        -html \
        -out $OUTPUT/blast_low_cov_${SAMPLE}_2_${REF}_vs_nt.html
done
