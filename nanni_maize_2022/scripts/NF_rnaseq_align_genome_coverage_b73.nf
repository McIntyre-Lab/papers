#!/usr/bin/env nextflow

/*
*    Nextflow script
*    to generate report use:  -with-report file-name
*    export _JAVA_OPTIONS=-Djava.io.tmpdir=/blue/mcintyre/share/maize_ainsworth/ROZ_NF
*/

// cat together techReps (different lanes) before aligning
// (1) align merged reads - SE alignment to B73 genome
// (2) align unmerged reads - PE alignemnt to B73 genome
// (3) combine SE and PE B73 genome alignments with samtools merge
// (4) coverage counts on fusions


println """\
    DIRECTORIES AND INPUT FILES
    input RG reads:       ${INPUT}
    input NG reads:       ${INPUTNG}
    output genome alns:   ${OUTG}
    scripts:              ${PY}
    design file:          ${DESIGN_FILE}
    """
    .stripIndent()


// split design file into 1 row chunks
// want to execute a task for each row
Channel
    .fromPath( DESIGN_FILE )
    .splitCsv( header: ['SAMPLE'] )
    .set { chunks_ch }


process aln_merged_SE_2_genome {

    input:
    val row from chunks_ch 

    output: 
    val row into genome_ch

    shell:
    '''
        
    module load bwa/0.7.17 gcc/9.3.0 python/2 samtools/1.10 bowtie/1.2.3

    mkdir -p !{OUTG}
    
    ROZ=!{PROJ}/roz_!{row.SAMPLE}
        mkdir -p ${ROZ}

    for i in bbmerge R1_unmerge R2_unmerge
    do
        ## cat lanes and remove empty reads
        cat !{INPUT}/!{row.SAMPLE}_*_${i}_min*_distinct.fq !{INPUTNG}/!{row.SAMPLE}_${i}_min*_distinct.fq > ${ROZ}/!{row.SAMPLE}_${i}.fq
    done

    echo "SE alignment"
    READ=${ROZ}/!{row.SAMPLE}_bbmerge.fq
    echo " read is:  ${READ}"
    echo " ref is:  !{GREF}/Zea_mays.B73_RefGen_v4.dna_sm.toplevel_BWA"

    bwa mem \
        -t 8 \
        -M !{GREF}/Zea_mays.B73_RefGen_v4.dna_sm.toplevel_BWA \
        ${READ} \
        >!{OUTG}/!{row.SAMPLE}_bbmerge_bwa_SE.sam

    echo " done with SE alignment to genome "

    '''
}


process aln_unmerged_PE_2_genome {

    input:
    val row from genome_ch

    output:
    val row into parse_ch

    shell:
    '''
    module load bwa/0.7.17 gcc/9.3.0 python/2 samtools/1.10 bowtie/1.2.3

    ROZ=!{PROJ}/roz_!{row.SAMPLE}

    READ1=${ROZ}/!{row.SAMPLE}_R1_unmerge.fq
    READ2=${ROZ}/!{row.SAMPLE}_R2_unmerge.fq

    echo "PE alignment"
    echo "READS are:  ${READ1} and ${READ2}"
    echo "reference is:  !{GREF}/Zea_mays.B73_RefGen_v4.dna_sm.toplevel_BWA"

    bwa mem \
        -t 8 \
        -M !{GREF}/Zea_mays.B73_RefGen_v4.dna_sm.toplevel_BWA \
        ${READ1} \
        ${READ2} \
        >!{OUTG}/!{row.SAMPLE}_unmerge_bwa_PE.sam

    echo " done with PE alignment to genome "

    '''
}

process parse_aln {

    input:
    val row from parse_ch

    output:
    val row into prep_ch

    shell:
    '''

    module load python/2

    mkdir -p !{PSAM}
    ROZ=!{PROJ}/roz_!{row.SAMPLE}
        mkdir -p ${ROZ}

    echo "parsing sam alignments! " 
    python /blue/mcintyre/share/python.git/BWASplitSAM_07mai.py \
        -fq1 ${ROZ}/!{row.SAMPLE}_bbmerge.fq \
        -s ${OUTG}/!{row.SAMPLE}_bbmerge_bwa_SE.sam \
        --outdir !{PSAM}
  
    python /blue/mcintyre/share/python.git/BWASplitSAM_07mai.py \
        -fq1 ${ROZ}/!{row.SAMPLE}_R1_unmerge.fq \
        -fq2 ${ROZ}/!{row.SAMPLE}_R2_unmerge.fq \
        -s !{OUTG}/!{row.SAMPLE}_unmerge_bwa_PE.sam \
        --outdir !{PSAM}

    echo "done parsing alignments"
    '''
}

process prep4counts {

    input:
    val row from prep_ch

    output:
    val row into cnts_ch

    shell:
    '''
    mkdir -p !{BAMFILES}

    module load bwa/0.7.17 gcc/9.3.0 python/2 samtools/1.10 bowtie/1.2.3
    
    echo "sam to bam, sort and index assembled and unassembles
    "
    REF1=!{GREF}/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa
    echo "samtools ref:  ${REF1} 
    "
    
    for i in bbmerge_bwa_SE unmerge_bwa_PE
    do

        echo "catting mapped and oposite"
        cat !{PSAM}/!{row.SAMPLE}_${i}_mapped.sam !{PSAM}/!{row.SAMPLE}_${i}_oposite.sam > !{PSAM}/!{row.SAMPLE}_${i}_uniq.sam
        samtools view -ut ${REF1}.fai -o !{BAMFILES}/!{row.SAMPLE}_${i}_uniq.bam !{PSAM}/!{row.SAMPLE}_${i}_uniq.sam

        samtools sort -T !{BAMFILES}/!{row.SAMPLE}_${i}_uniq.tmp.sorted \
            -o !{BAMFILES}/!{row.SAMPLE}_${i}_uniq.sorted.bam \
            -O bam !{BAMFILES}/!{row.SAMPLE}_${i}_uniq.bam
        samtools index !{BAMFILES}/!{row.SAMPLE}_${i}_uniq.sorted.bam
    done


    echo "combine SE and PE bam files
    "
    samtools merge -f \
        !{BAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.bam \
        !{BAMFILES}/!{row.SAMPLE}_unmerge_bwa_PE_uniq.bam \
        !{BAMFILES}/!{row.SAMPLE}_bbmerge_bwa_SE_uniq.bam

    samtools sort -T !{BAMFILES}/!{row.SAMPLE}_combined_bwa.tmp.sorted \
        -o !{BAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.bam \
        -O bam !{BAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.bam
    samtools index !{BAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.bam

    '''
}

 
process cvr {

    input:
    val row from cnts_ch

    output:
    val row into clean_ch

    shell:
    '''
    mkdir -p !{CVR}
    mkdir -p !{LOGC}
    mkdir -p !{SAMFILES}
    mkdir -p !{PILES}

    module load bwa/0.7.7 gcc/5.2.0 python/2.7.10 samtools

    FASTA=!{GREF}/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa
    ROZ=!{PROJ}/roz_!{row.SAMPLE}

    ## create mpileup from bam files
    samtools index !{BAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.bam
    samtools mpileup -d 1000000 -f ${FASTA} !{BAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.bam \
        > !{PILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.mpileup


    ## convert sorted bam to sorted sam
    samtools view -h -o !{SAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.sam !{BAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.bam

    ## get coverage counts for alignments using rpkm_calculate
    python !{PROJ}/scripts/rpkm_calculate.py \
        -m !{PILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.mpileup \
        -s !{SAMFILES}/!{row.SAMPLE}_combined_bwa_uniq.sorted.sam \
        -b !{BED} \
        -g !{LOGC}/cvrg_cnts_!{row.SAMPLE}.logfile \
        -o ${ROZ}/cvrg_cnts_!{row.SAMPLE}.csv \
        -c

    ## add column with sampleID to cov counts
    awk -F "," -v id="!{row.SAMPLE}" -v OFS="," '{print id, $0}' ${ROZ}/cvrg_cnts_!{row.SAMPLE}.csv > !{CVR}/cvrg_cnts_!{row.SAMPLE}.csv

    '''    
}


process cleanUp {

    input:
    val row from clean_ch.collect()

    shell:
    '''
#    rm -r ${ROZ}
#    rm -r !{OUTG}
#    rm !{PSAM}/*.sam
#    rm !{PSAM}/*.fq
#    rm -r !{PSAM}/mpileups
#    rm -r !{PSAM}/sam_files_uniq

    '''
}

