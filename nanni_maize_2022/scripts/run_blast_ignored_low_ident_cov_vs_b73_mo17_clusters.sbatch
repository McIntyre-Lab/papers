#!/bin/bash

#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=blast
#SBATCH --output=/blue/mcintyre/share/maize_ainsworth/scripts/SLURM_LOGS/blast_2_clusters_%A_%a.out
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH --cpus-per-task=1
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre
#SBATCH --array=4,6-9

## Non-B73/Mo17 samples are arrays 3,4,6-9

### blast PB transcripts low identity and low coverage to b73 and mo17 clusters

module load ncbi_blast python seqtk

### Set Directories
PROJ=/blue/mcintyre/share/maize_ainsworth
OUTPUT=$PROJ/compare_b73_2_mo17/blast_low_identity_PBreads
    mkdir -p ${OUTPUT}

## Design file
DESIGN_FILE=$PROJ/design_files/df_maize_test_PacBio_fullpath_noHeader_allSamples.csv
DESIGN=$(cat $DESIGN_FILE | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)
IFS=',' read -ra ARRAY <<< "$DESIGN"

GENO=${ARRAY[1]}
TRT=${ARRAY[2]}    
ID=${ARRAY[3]}

SAMPLE=${ID}_${GENO}_${TRT}

## Set up temp directory
ROZ=${OUTPUT}/roz_${SAMPLE}
    mkdir -p ${ROZ}

for REF in b73; do
#for REF in b73 mo17_yan mo17_cau; do
    
    for TYPE in ident cov; do
        ## Set output file
        if [[ ${TYPE} == "ident" ]]; then
            OUTPUT=$PROJ/compare_b73_2_mo17/blast_low_identity_PBreads
        else
            OUTPUT=$PROJ/compare_b73_2_mo17/blast_low_coverage_PBreads
        fi
        mkdir -p ${OUTPUT}

        FA=$PROJ/compare_b73_2_mo17/uncollapsed_${REF}_IsoSeq_reads/${SAMPLE}_unCollapsed_low_${TYPE}_${REF}_ref.fa

        for COMPARE in b73 mo17; do
            ### query = low identity fasta 
            ### subject = nr database

            BLASTDB=$PROJ/compare_b73_2_mo17/blast_db_collapsed_clusters/${COMPARE}_collapsed.fa.BLAST

            date
            echo ${FA}
            echo ${SAMPLE}
            echo ${BLASTDB}

            QUERY=${FA}

            ## create header for blast tsv output
            echo "qseqid	stitle	qstart	qend	qlen	sstart	send	slen	length	bitscore	evalue	mismatch	gaps" > $ROZ/header.tsv

            ## blast and output csv file
            blastn \
                -db ${BLASTDB} \
                -query $QUERY \
                -outfmt "6 qseqid stitle qstart qend qlen sstart send slen length bitscore evalue mismatch gaps" \
                -out ${ROZ}/blast_low_${TYPE}_${SAMPLE}_2_${REF}_vs_${COMPARE}_collapsed_clusters.tsv

            ## add header line
            cat ${ROZ}/header.tsv ${ROZ}/blast_low_${TYPE}_${SAMPLE}_2_${REF}_vs_${COMPARE}_collapsed_clusters.tsv > ${OUTPUT}/blast_low_${TYPE}_${SAMPLE}_2_${REF}_vs_${COMPARE}_collapsed_clusters.tsv

            ### blast and output html file
            blastn \
                -db ${BLASTDB} \
                -query ${QUERY} \
                -html \
                -out ${OUTPUT}/blast_low_${TYPE}_${SAMPLE}_2_${REF}_vs_${COMPARE}_collapsed_clusters.html
        done
    done
done

rm -r ${ROZ}
