#!/bin/sh

### EVENT ANALYSIS EVENT SUMMARY AND TRANSCRIPT ELIMINATION WORKFLOW

# Set paths and references
PROJ=/mnt/store/event_sandbox/python_testing

# Design file
DESIGN=$PROJ/input_files/NPC_OPC_design_file.tsv

# Input counts
FRAGCOUNTS=${PROJ}/coverage_counts_fragments
FUSCOUNTS=${PROJ}/coverage_counts_fusions
INTCOUNTS=${PROJ}/coverage_counts_introns
JUNCCOUNTS=${PROJ}/coverage_counts_splicing

# Output path
OUTDIR=$PROJ/output_files
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

# Set detection parameters
MINAPN=0
MINDTCT=0.5

# Set grouping variable
GROUPVAR=cell_type

### (1) Convert single-file coverage counts into a single wide dataset
echo "Creating wide coverage datasets"
# Fragments 
python $PROJ/programs/import_counts_and_convert_to_wide.py --input-directory ${FRAGCOUNTS} \
                                            --output-wide ${OUTDIR}/NPC_OPC_fragment_counts_wide.tsv

# Fusions (exonic regions) 
python $PROJ/programs/import_counts_and_convert_to_wide.py --input-directory ${FUSCOUNTS} \
                                            --output-wide ${OUTDIR}/NPC_OPC_fusion_counts_wide.tsv

# Introns
python $PROJ/programs/import_counts_and_convert_to_wide.py --input-directory ${INTCOUNTS} \
                                            --output-wide ${OUTDIR}/NPC_OPC_intron_counts_wide.tsv

# Junctions 
python $PROJ/programs/import_counts_and_convert_to_wide.py --input-directory ${JUNCCOUNTS} \
                                            --output-wide ${OUTDIR}/NPC_OPC_junction_counts_wide.tsv


### (2) Flag events on/off
echo "Identifying detected events"

# Fragments 
python $PROJ/programs/flag_event_detection.py --input-data ${OUTDIR}/NPC_OPC_fragment_counts_wide.tsv \
                                              --design-file ${DESIGN} \
                                              --group-variable ${GROUPVAR} \
                                              --minimum-abundance ${MINAPN} \
                                              --minimum-proportion ${MINDTCT} \
                                              --output-flags ${OUTDIR}/NPC_OPC_fragment_flagged_apn0.tsv

# Fusions (exonic regions) 
python $PROJ/programs/flag_event_detection.py --input-data ${OUTDIR}/NPC_OPC_fusion_counts_wide.tsv \
                                              --design-file ${DESIGN} \
                                              --group-variable ${GROUPVAR} \
                                              --minimum-abundance ${MINAPN} \
                                              --minimum-proportion ${MINDTCT} \
                                              --output-flags ${OUTDIR}/NPC_OPC_fusion_flagged_apn0.tsv

# Introns
python $PROJ/programs/flag_event_detection.py --input-data ${OUTDIR}/NPC_OPC_intron_counts_wide.tsv \
                                              --design-file ${DESIGN} \
                                              --group-variable ${GROUPVAR} \
                                              --minimum-abundance ${MINAPN} \
                                              --minimum-proportion ${MINDTCT} \
                                              --output-flags ${OUTDIR}/NPC_OPC_intron_flagged_apn0.tsv

# Junctions 
python $PROJ/programs/flag_event_detection.py --input-data ${OUTDIR}/NPC_OPC_junction_counts_wide.tsv \
                                              --design-file ${DESIGN} \
                                              --group-variable ${GROUPVAR} \
                                              --minimum-abundance ${MINAPN} \
                                              --minimum-proportion ${MINDTCT} \
                                              --output-flags ${OUTDIR}/NPC_OPC_junction_flagged_apn0.tsv

### (3) Create event-level summaries
echo "Creating event-level summaries"
# Fragments 
python $PROJ/programs/create_event_summaries.py --input-data ${OUTDIR}/NPC_OPC_fragment_counts_wide.tsv \
                                                --design-file ${DESIGN} \
                                                --event-length 10 \
                                                --group-variable ${GROUPVAR} \
                                                --annotation-file ${OUTDIR}/mm10_exon_fragment_annotations.csv \
                                                --detection-flags-file ${OUTDIR}/NPC_OPC_fragment_flagged_apn0.tsv \
                                                --output-summary ${OUTDIR}/NPC_OPC_summary_by_exon_fragment.tsv

# Fusions (exonic regions) 
python $PROJ/programs/create_event_summaries.py --input-data ${OUTDIR}/NPC_OPC_fusion_counts_wide.tsv \
                                                --design-file ${DESIGN} \
                                                --event-length 10 \
                                                --group-variable ${GROUPVAR} \
                                                --annotation-file ${OUTDIR}/mm10_fusion_annotations.csv \
                                                --detection-flags-file ${OUTDIR}/NPC_OPC_fusion_flagged_apn0.tsv \
                                                --output-summary ${OUTDIR}/NPC_OPC_summary_by_fusion.tsv

# Introns
python $PROJ/programs/create_event_summaries.py --input-data ${OUTDIR}/NPC_OPC_intron_counts_wide.tsv \
                                                --design-file ${DESIGN} \
                                                --group-variable ${GROUPVAR} \
                                                --event-length 10 \
                                                --annotation-file ${OUTDIR}/mm10_introns_from_fusions.csv \
                                                --detection-flags-file ${OUTDIR}/NPC_OPC_intron_flagged_apn0.tsv \
                                                --output-summary ${OUTDIR}/NPC_OPC_summary_by_intron.tsv

# Junctions 
python $PROJ/programs/create_event_summaries.py --input-data ${OUTDIR}/NPC_OPC_junction_counts_wide.tsv \
                                                --design-file ${DESIGN} \
                                                --event-length 56 \
                                                --group-variable ${GROUPVAR} \
                                                --junction-sequence-index ${OUTDIR}/mm10_junction_to_sequence_index.csv \
                                                --annotation-file ${OUTDIR}/mm10_junction_annotations.csv \
                                                --detection-flags-file ${OUTDIR}/NPC_OPC_junction_flagged_apn0.tsv \
                                                --output-summary ${OUTDIR}/NPC_OPC_summary_by_junction.tsv

### (4) Identify expressed genes in each group
echo "Identifying expressed genes"
python $PROJ/programs/identify_expressed_genes.py --input-exonic-region-summary ${OUTDIR}/NPC_OPC_summary_by_fusion.tsv \
                                                  --design-file ${DESIGN} \
                                                  --group-variable ${GROUPVAR} \
                                                  --output-gene-summary ${OUTDIR}/NPC_OPC_summary_of_expressed_genes.tsv

### (5) Create intron-border junction summaries
echo "Creating intron-to-border junction summaries"
python $PROJ/programs/create_intron_border_summary.py \
       --input-exonic-region-summary ${OUTDIR}/NPC_OPC_summary_by_fusion.tsv \
       --input-intron-summary ${OUTDIR}/NPC_OPC_summary_by_intron.tsv \
       --input-junction-summary ${OUTDIR}/NPC_OPC_summary_by_junction.tsv \
       --intron-border-junction-index ${OUTDIR}/mm10_intron2border_junction_index.csv \
       --design-file ${DESIGN} \
       --group-variable ${GROUPVAR} \
       --minimum-donor-mean 5 \
       --minimum-acceptor-mean 5 \
       --minimum-intron-mean 5 \
       --output-intron-border-summary ${OUTDIR}/NPC_OPC_summary_of_intron_borderjunctions.tsv

outSummaryoutSummary
### (6) Summarize transcripts
echo "Summarizing transcripts"
python $PROJ/programs/summarize_transcripts_by_group.py --input-gene-summary ${OUTDIR}/NPC_OPC_summary_of_expressed_genes.tsv \
--input-exon-fragment-summary ${OUTDIR}/NPC_OPC_summary_by_exon_fragment.tsv \
--input-junction-summary ${OUTDIR}/NPC_OPC_summary_by_junction.tsv \
--input-event-to-transcript-index ${OUTDIR}/mm10_event2transcript2gene_index.csv \
--design-file ${DESIGN} \
--group-variable ${GROUPVAR} \
--output-transcript-summary ${OUTDIR}/NPC_OPC_summary_of_transcripts_exp_genes.tsv

