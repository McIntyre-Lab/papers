#!/bin/sh

### EVENT ANALYSIS ANNOTATION CREATION WORKFLOW

# Set paths and references
PROJ=/mnt/store/event_sandbox/python_testing
GFF=$PROJ/input_files/aconesa_refseq_gff3_v4.gff
FASTA=$PROJ/input_files/mm10_for_bedtools_v2.fa

# Output path
OUTDIR=$PROJ/output_files
if [ ! -e $OUTDIR ]; then mkdir -p $OUTDIR; fi

### (1) EXTRACT EVENTS
# Extract exons
echo "Extracting exons"
#python2.7 $PROJ/programs/extractExons.py --gff ${GFF} --obed $OUTDIR/mm10_exons.bed --otable $OUTDIR/mm10_exons.csv

# Extract fusions (exonic regions)
echo "Extracting fusions"
#python2.7 $PROJ/programs/buildFusions.py --gff ${GFF} --obed $OUTDIR/mm10_fusions.bed --otable $OUTDIR/mm10_fusions.tsv

# Extract exon fragments
echo "Extracting exon fragments"
#python2.7 $PROJ/programs/extractExonFragments.py --gff ${GFF} --output $OUTDIR/mm10_exon_fragments.csv

# Build unambiguous introns from fusions
echo "Building unambiguous introns"
#python $PROJ/programs/build_unambiguous_introns_from_fusions.py --input-fusion-file $OUTDIR/mm10_fusions.tsv \
#                                                                --input-fusion-bed $OUTDIR/mm10_fusions.bed \
#                                                                --outCSV $OUTDIR/mm10_introns_from_fusions.csv \
#                                                                --outBED $OUTDIR/mm10_introns_from_fusions.bed

# Extract logical junctions
echo "Extracting logical junctions"
#python2.7 $PROJ/programs/extractJunctions.py --input ${GFF}.db --output $OUTDIR/mm10_logical_junctions.bed --size 40

# Extract junctions annotated to transcripts
echo "Extracting annotated junctions"
#python2.7 $PROJ/programs/extractTranscriptJunctions.py --gff ${GFF} --output $OUTDIR/mm10_annotated_junctions.csv

# Extract exon-skipping junction annotations
echo "Extracting exon-skipping annotations"
#python2.7 $PROJ/programs/extractSkippedExons.py --gff ${GFF} --otable $OUTDIR/mm10_exonskipping_junctions.csv \
#                                                --olist $OUTDIR/mm10_skipped_exons_list.csv

# Build donor-site exon-intron border junctions
echo "Building donor-site exon-intron border junctions"
#python2.7 $PROJ/programs/build_donor_border_junctions.py --gff ${GFF} --obed $OUTDIR/mm10_donor_border_junctions.bed \
#                                       --otable $OUTDIR/mm10_donor_border_junctions.csv --size 40

# Build acceptor-site exon-intron border junctions
echo "Building acceptor-site exon-intron border junctions"
#python2.7 $PROJ/programs/build_acceptor_border_junctions.py --gff ${GFF} \
#                                       --obed $OUTDIR/mm10_acceptor_border_junctions.bed \
#                                       --otable $OUTDIR/mm10_acceptor_border_junctions.csv --size 40

### (2) EXON ANNOTATIONS
# Format exon annotations
echo "Format exon annotations"
python $PROJ/programs/import_and_format_exons.py --input $OUTDIR/mm10_exons.csv --output $OUTDIR/mm10_exon_annotations.csv

# Format fusion (exonic region) annotations
echo "Format fusion annotations"
python $PROJ/programs/import_and_format_fusions.py --input-fusion-file $OUTDIR/mm10_fusions.tsv \
                                    --input-fusion-bed $OUTDIR/mm10_fusions.bed \
                                    --input-exon-file $OUTDIR/mm10_exons.csv \
                                    --outCSV $OUTDIR/mm10_fusion_annotations.csv \
                                    --outBED $OUTDIR/mm10_fusions_coverage.bed

# Format exon fragment annotations
echo "Format exon fragment annotations"
python $PROJ/programs/import_and_format_fragments.py --input-fragment-file $OUTDIR/mm10_exon_fragments.csv \
                                      --input-fusion-file $OUTDIR/mm10_fusions.tsv \
                                      --input-fusion-bed $OUTDIR/mm10_fusions.bed \
                                      --input-exon-file $OUTDIR/mm10_exons.csv \
                                      --outCSV $OUTDIR/mm10_exon_fragment_annotations.csv \
                                      --outBED $OUTDIR/mm10_exon_fragments_coverage.bed

### (3) JUNCTION ANNOTATIONS
# Format logical junctions
echo "Format logical junctions"
python $PROJ/programs/import_and_format_junctions.py --bed $OUTDIR/mm10_logical_junctions.bed \
                                      --output $OUTDIR/mm10_logical_junctions_formatted.csv

# Flag annotated junctions
echo "Flag annotated junctions"
python $PROJ/programs/import_and_flag_transcript_junctions.py --input-junctions $OUTDIR/mm10_logical_junctions_formatted.csv \
                                               --input-annotated-junctions $OUTDIR/mm10_annotated_junctions.csv \
                                               --output $OUTDIR/mm10_logical_junctions_flag_annotated.csv

# Flag exon-skipping junctions
echo "Add exon-skipping annotations"
python $PROJ/programs/import_exon_skipping_annotations.py --input-junctions $OUTDIR/mm10_logical_junctions_flag_annotated.csv \
                                           --input-exonskip-annot $OUTDIR/mm10_exonskipping_junctions.csv \
                                           --output $OUTDIR/mm10_logical_junctions_flag_exonskip.csv

# Append border junctions
echo "Append border junctions"
python $PROJ/programs/append_border_junctions.py --input-junctions $OUTDIR/mm10_logical_junctions_flag_exonskip.csv \
                                  --input-donor-border-junctions $OUTDIR/mm10_donor_border_junctions.bed \
                                  --input-acceptor-border-junctions $OUTDIR/mm10_acceptor_border_junctions.bed \
                                  --junction-size 40 \
                                  --output $OUTDIR/mm10_logical_junctions_and_border_junctions.csv

# Add exon info to junctions
echo "Add exon information to junctions"
python $PROJ/programs/add_exon_info_to_junctions.py --input-junction-file $OUTDIR/mm10_logical_junctions_and_border_junctions.csv \
                                     --input-exon-file $OUTDIR/mm10_exon_annotations.csv \
                                     --output-junction-info $OUTDIR/mm10_junctions_w_exon_info.csv

# Collapse junctions with identical coordinates
echo "Collapse junctions with identical coordinates"
python $PROJ/programs/collapse_duplicate_junctions.py --input-junction-file $OUTDIR/mm10_junctions_w_exon_info.csv \
                                       --input-exon-annotation $OUTDIR/mm10_exons.csv \
                                       --output-collapsed-junctions $OUTDIR/mm10_junctions_full_annotation.csv

# Extract junction sequences
echo "Extracting junction sequences"
python $PROJ/programs/extract_junction_sequence.py --input-junction-file $OUTDIR/mm10_junctions_full_annotation.csv \
                                    --input-fasta-file ${FASTA} \
                                    --output-junction-annotation $OUTDIR/mm10_junction_annotations.csv \
                                    --output-junction-to-seq-index $OUTDIR/mm10_junction_to_sequence_index.csv \
                                    --output-junction-sequences $OUTDIR/mm10_junctions.fa \
                                    --output-coverage-bed $OUTDIR/mm10_junctions_coverage.bed 

### (4) CREATE EVENT INDICES
# Make event-to-transcript-to-gene index
echo "Creating event-to-transcript-gene index"
python $PROJ/programs/build_Event2Transcript_index.py -e $OUTDIR/mm10_exons.csv \
                                       -f $OUTDIR/mm10_exon_fragment_annotations.csv \
                                       -j $OUTDIR/mm10_junction_annotations.csv \
                                       -o $OUTDIR/mm10_event2transcript2gene_index.csv

# Make intron-to-border junction index
echo "Creating intron-to-border-junction index"
python $PROJ/programs/build_intron2border_junction_index.py --intron-annotation-file $OUTDIR/mm10_introns_from_fusions.csv \
                                             --junction-annotation-file $OUTDIR/mm10_junction_annotations.csv \
                                             --output-intron-index $OUTDIR/mm10_intron2border_junction_index.csv

s
