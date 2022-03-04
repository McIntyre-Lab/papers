#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=blast
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/blast_%A_%a.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
##SBATCH --array=1-11
#SBATCH --array=10

### blast all CCS that do not map to b73, mo17_yan, or mo17_cau to ncbi nr database

module load ncbi_blast/2.10.1 python/3.6 samtools/1.10

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
IND=$PROJ/check_isoseq3_lost_genes/minimap2
OUTPUT=$PROJ/blast_unmapped_ccs
    mkdir -p $OUTPUT

ROZ=$IND/rozBlast
    mkdir -p $ROZ


# Get info from design file

DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_path_noHeader_all.csv
DESIGN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $DESIGN_FILE)
IFS=',' read -ra ARRAY <<< "$DESIGN"

PACBIO_NUM=${ARRAY[0]}
GENO=${ARRAY[1]}
TRT=${ARRAY[2]}
ID=${ARRAY[3]}
SAMPLE_PATH=${ARRAY[4]}

SAMPLE_NAME=${ID}_${GENO}_${TRT}


for REF in b73 mo17_yan mo17_cau
do

    ### subject = nr database

    ## Select unmapped reads from SAM (where sam flag is 4) and convert to FASTA
    SAM=${IND}/${SAMPLE_NAME}_unfiltered_minimap2_${REF}.sam
    FA=${IND}/${SAMPLE_NAME}_unmapped_ccs_2_${REF}.fa

    samtools view -f 4 ${SAM} | samtools fasta > ${FA}

    echo ${SAMPLE_NAME}
    echo ${FA}

    QUERY=${FA}
        
    ## create header for blast tsv output
    echo "qseqid	stitle	qstart	qend	qlen	sstart	send	slen	length	bitscore	evalue	mismatch	gaps" > $ROZ/header.tsv

    ## blast and output csv file
    blastn \
        -db nt \
        -query $QUERY \
        -outfmt "6 qseqid stitle qstart qend qlen sstart send slen length bitscore evalue mismatch gaps" \
        -out $ROZ/${SAMPLE_NAME}_blast_unmapped_vs_nt.tsv

    ## add header line
    cat $ROZ/header.tsv $ROZ/${SAMPLE_NAME}_blast_unmapped_vs_nt.tsv \
        > $OUTPUT/${SAMPLE_NAME}_blast_unmapped_ccs_vs_nt.tsv

    ### blast and output html file
    blastn \
        -db nt \
        -query $QUERY \
        -html \
        -out $OUTPUT/${SAMPLE_NAME}_blast_unmapped_ccs_vs_nt.html
done
