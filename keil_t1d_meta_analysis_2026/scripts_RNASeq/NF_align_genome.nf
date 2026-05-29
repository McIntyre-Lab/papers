#!/usr/bin/env nextflow

/*
*    Nextflow script for RNAseq alignment check 
*    
*    align dataPrepped reads to genome
*
*    to generate report use:  -with-report file-name
*    export _JAVA_OPTIONS=-Djava.io.tmpdir=/blue/concannon/share/mcintyre/ROZ_NF
*
*/


// Directories
PROJ    =   "/blue/concannon/share/T1DGC_R01_controls_RNAseq"
SCRIPTS =   "${PROJ}/scripts"
REF     =   "/blue/concannon/share/jnewman/references/hg38_ensembl"
READS   =   "${PROJ}/dataPrep_noDup_uniq_reads"
ALN     =   "${PROJ}/alignments_hg38"
DF      =   "${PROJ}/design_files"

// Files
DESIGN_FILE =   "${DF}/design_all_sampleID_techRep.csv"
// sampleID,dir,techRep
// 40359602-CD19,1910UNHS-0108,1
// 40359602-CD19,2005UAHS-0081,2

println """\
    DIRECTORIES AND INPUT FILES
    input reads:        ${READS}
    output alignments:  ${ALN}
    scripts:            ${SCRIPTS}
    design file:        ${DESIGN_FILE}
    """
    .stripIndent()


// split design file into 1 row chunks 
// want to execute a task for each row
Channel
    .fromPath( DESIGN_FILE )
    .splitCsv( header: ['DIR','SAMPLEID','TR'], skip: 1 )
    .set { chunks_ch }


process BWA {

    input:
    val row from chunks_ch

    shell:

    '''
    mkdir -p !{ALN}

    module load bwa/0.7.7 samtools/1.12

    ### BWA alignments
    bwa mem -t 8 \
        -M \
	!{REF}/hg38.nochr.fa \
        !{READS}/!{row.SAMPLEID}_!{row.TR}_noDup_uniq.fastq \
        > !{ALN}/!{row.SAMPLEID}_!{row.TR}.sam

    ## convert to bam and remove sam
    samtools view -h -S -b !{ALN}/!{row.SAMPLEID}_!{row.TR}.sam > !{ALN}/!{row.SAMPLEID}_!{row.TR}.bam

    rm !{ALN}/!{row.SAMPLEID}_!{row.TR}.sam

    '''
}

