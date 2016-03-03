/* Makefile for creating the splicing catalogue */

/* Set path to SAS programs */

/* Import the generated output of all possible, logical junctions within a gene */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/import_junctions_and_format.sas"

/* Import exon information */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/import_exons_and_format.sas"

/* Import junctions generated from only transcripts */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/import_junctions_from_transcripts.sas"

/* Import annotations of skipped exons */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/import_skipped_exon_annotations.sas"

/* Import generated intron retention events */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/import_intron_retention_events.sas"

/* Annotated junctions if they exist in transcripts or not */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/flag_annotated_junctions.sas"

/* Annotate junctions if exon-skipping */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/flag_exonskip_junctions.sas"

/* Stack junctions and intron retention events */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/cat_junctions_and_ir.sas"

/* Add exon information to splicing events and annotate alternative donors and acceptors */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/add_exon_info.sas"

/* Merge duplicate junctions based on coordinate */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/collapse_duplicates.sas"

/* Format annotations of splicing events */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/format_splicing_annotations.sas"

/* Generate BED file of splicing events */
    %include "/home/jrbnewman/McLab/junction_annotations/sas_programs/hg19/catalogue2BED.sas"


