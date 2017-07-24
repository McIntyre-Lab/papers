#!/bin/sh
#########################################################################################
# 
# Title: EVENT ANALYSIS ANNOTATIONS MAKEFILE : PART 1, BUILDING ANNOTATIONS
# Author: Jeremy R. B. Newman
# 
# Description: This script runs through all the associated Python and SAS programs to
# create the annotation files needed by Event Analysis to identify the probable
# transcriptome. It takes a GTF/GFF file as input. The output files are a series of CSVs
# that are then imported into SAS to build the annotation datasets.
# 
# There are three main outputs:
# (1) Fusion (exonic regions) -- these are regions of overlapping exons derived
#     from the reference annotations
# (2) Exon fragments -- this takes the regions of overlapping exons are divides them
#     based on their transcriptional identity
# (3) All possible, logical junctions that can be derived from the exons of a gene
# 
# This script requires the following software to be installed:
#	Python 2.7.6 or above (test on Py3.5), with the following packages:
#               gffutils 0.8.5 or above
#               itertools
#		mclab libraries (github link)
#	BEDtools
#       SAS 9.3 or above
#
# You will need to modify some parameters within the script, such as GTF/GFF location,
# read length and project folder (see SECTION 1 below).
#
# Usage: sh ./run_event_analysis_build_annotations.sh
#
#########################################################################################

#########################################################################################
#
#  NOTE ON GFF3 FORMATTING:
#
#  The scripts used in this pipeline were developed based on the formatting of the Flybase
#  D.mel release 5/6 GFF3 format, and requires the following GFF entries:
#      gene
#      transcript, mRNA, and/or ncRNA
#      exon
#      chromosome or chromosome_arm
#
#  Each entry should have an ID or Name, and Parent (if applicable) attributes. For example:
#      For "gene" entries:                      ID=UBASH3A;Name=UBASH3A
#      For "transcript" entries:                ID=UBASH3A.gAug10;Name=UBASH3A.gAug10;Parent=UBASH3A
#      For "exon" entries                       ID=UBASH3A:20;Name=UBASH3A:20;Parent=UBASH3A.aAug10,UBASH3A.bAug10
#                                               (Unique features should have a single entry, even if
#                                               it has multiple parent features)
#
#  Even if your GFF looks like this, we highly recommend the convertGTF2GFF3.py script be
#  run to ensure that the GFF entries are correct.
#
#  NOTE: Your input GTF/GFF _MUST_ be sorted by chromsome, start position and stop position!
#        Otherwise you will likely run into errors while creating annotations.
#
#        If you need to sort your GTF/GFF file, run the following in terminal:
#
#              sort -k1n -k4n -k5n $GFF > $GFF.sorted
#
#         Then use this sorted GFF as the input GTF/GFF as input for the scripts below
#
###########################################################################################

###########################################################################################
#
# SECTION 1: VARIABLES TO BE SET BY USER
#
###########################################################################################

echo "Setting project folder"
## Set project folder path. This will be where you want to store the annotation files
PROJ=/mnt/store/event_sandbox/makefile_testing
   if [ ! -e $PROJ ]; then mkdir -p $PROJ; fi

echo "Setting SAS library folder"
## SAS library folder will be the location where the SAS datasets are stored
SASLIB=$PROJ/sas_data
   if [ ! -e $SASLIB ]; then mkdir -p $SASLIB; fi

### Set prefix for output files. This should be short but descriptive of the genome build
### you are using, e.g. fb611 for Flybase 6.11
PREFIX=fb611

echo "Setting path to scripts and SAS programs"
## Set the location of scripts and SAS programs necessary to run the data
SCRIPTS=$PROJ/scripts
SASPROG=$PROJ/sas_programs
OUTDIR=$PROJ/generated_files

   if [ ! -e $SCRIPTS ]; then mkdir -p $PROJ; fi
   if [ ! -e $SASPROG ]; then mkdir -p $PROJ; fi
   if [ ! -e $OUTDIR ]; then mkdir -p $PROJ; fi

echo "Setting SAS temp directory"
## SAS temp/scratch directory location
SASWORK=/mnt/SAS_WRK2/sas_work
   if [ ! -e $SASWORK ]; then mkdir -p $SASWORK; fi

echo "Setting SAS log directory"
## Set the location of the catalog scripts
SASLOGS=$SASPROG/logs
   if [ ! -e $SASLOGS ]; then mkdir -p $SASLOGS; fi

echo "Setting GTF/GFF file"
## Set the path to your input GTF/GFF3 file
GFF=/home/jrbnewman/McLab/useful_dmel_data/flybase611/gff/dmel_r6.11_no_analysis_converted_v2.gff

echo "Setting junction size"
## Set read length. This will be the length of your sequencing reads
## The maximum size of each splicing event will be 24bp larger than the length of the read
## The EVENT_SIZE variable is half of this maximum length 
READ_LENGTH=50
JUNC_SIZE=$(expr ${READ_LENGTH} / 2 + 12)

echo "Read length is ${READ_LENGTH} bp"
echo "Setting junction size to $(expr ${JUNC_SIZE} \* 2) bp"

echo "Setting location of FASTA file"
## Set the path to the genome FASTA file
FASTA=/home/jrbnewman/McLab/useful_dmel_data/flybase611/flybase_files/dmel-all-chromosome-r6.11.fasta

###########################################################################################
#
# SECTION 2: ANNOTATION FILE GENERATION
#
###########################################################################################

# Checking if a gff.db file exists, and if not then create one
echo "Checking if user-supplied GTF/GFF file has a pre-generated database file"
    if [ ! -e ${GFF}.db ]
    then
    echo "No database file found! Generating..."
    python $SCRIPTS/make_gff_db.py --gff $GFF
    echo "Database generated!"
    else
    echo "Database file found! Skipping generation."
    fi

# Convert user-supplied GTF/GFF file into GFF3 format
echo "Converting user-supplied GTF/GFF to GFF3 format"
GFFOUT=${OUTDIR}/gff3_sorted_for_annot.gff3
python $SCRIPTS/convertGTF2GFF3.py --input $GFF --output ${GFFOUT}

python $SCRIPTS/make_gff_db.py --gff ${GFFOUT}

# Build fusions (exonic regions)
echo "Creating annotation CSVs for fusions (exonic regions)"
python $SCRIPTS/buildFusions.py --gff $GFFOUT \
                                --obed ${OUTPUT}/${PREFIX}_fusions_si.bed \
                                --otable ${OUTPUT}/${PREFIX}_fusions_si.tsv

# Build exon fragments
echo "Creating annotation CSVs for exon fragments"
python $SCRIPTS/exon_overlapping_pieces.py --gff $GFFOUT \
                                --obed ${OUTPUT}/${PREFIX}_exon_fragments.csv 

# Creating all possible, logical junctions within a gene, from combinations of existing exons
echo "Generating all possible, logical junctions from gene models..."
python $SCRIPTS/extractJunctions.py --input ${GFFOUT}.db \
                                    --output ${OUTPUT}/${PREFIX}_logical_junctions.bed \
                                    --size ${EVENT_SIZE}
echo "Done!"

# Making a list of all exons
echo "Generating list of all exons..."
python $SCRIPTS/extractExons.py --gff $GFFOUT \
                                --obed ${OUTPUT}/${PREFIX}_exons.bed \
                                --otable ${OUTPUT}/${PREFIX}_exons.csv
echo "Done!"

# Making a list of all junctions from transcripts
echo "Generating list of all junctions within transcripts..."
python $SCRIPTS/extractTranscriptJunctions.py --gff $GFFOUT \
                                              --output ${OUTPUT}/${PREFIX}_transcript_junctions.csv
echo "Done!"

# Making a list of all possible exon-skipping junctions
echo "Generating list of exon-skipping junctions..."
python $SCRIPTS/extractSkippedExons.py --gff $GFFOUT \
                                       --otable ${OUTPUT}/${PREFIX}_exon_skipping_junctions.csv \
                                       --olist ${OUTPUT}/${PREFIX}_skipped_exons_list.csv
echo "Done!"

# Making a list of all border junctions (exon-intron boundaries)
echo "Generating list of exon-intron border junctions..."
python $SCRIPTS/buildRetainedIntrons.py --gff $GFFOUT \
                                        --obed ${OUTPUT}/${PREFIX}_intron_retention.bed \
                                        --otable ${OUTPUT}/${PREFIX}_intron_retention.csv \
                                        --size ${EVENT_SIZE}
echo "Done!"

###########################################################################################
#
# SECTION 3: ANNOTATION FORMATTING AND PREPARATION
#
###########################################################################################


# Import exon, fusion/exonic region and fragment annotations
echo "Importing exon, exonic region and exon fragment annotations"
EXONINFO=${OUTPUT}/${PREFIX}_exons.csv
FUSBED=${OUTPUT}/${PREFIX}_fusions_si.bed
FUSINFO=${OUTPUT}/${PREFIX}_fusions_si.tsv
FRAGINFO=${OUTPUT}/${PREFIX}_exon_fragments.csv
EXONLEN=$(wc -l ${EXONINFO})
FUSBEDLEN=$(wc -l ${FUSBED})
FUSINFOLEN=$(wc -l ${FUSINFO})
FRAGLEN=$(wc -l ${FRAGINFO})

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/import_exon_fusion_fragment_annotations.sas \
    -work $SASWORK \
    -set libpath ${LIBPATH} \
    -set exonInfo ${EXONINFO} \
    -set fusBED ${FUSBED} \
    -set fusInfo ${FUSINFO} \
    -set fragInfo ${FRAGINFO} \
    -set prefix ${PREFIX}
    -set fusBEDLen ${FUSBEDLEN} \
    -set fusInfoLen ${FUSINFOLEN} \
    -set fragInfoLen ${FRAGLEN} \
    -set exonLen ${EXONLEN}

echo "Formatting exon annotations"
sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/format_exon_annotations.sas \
    -work $SASWORK \
    -set libpath ${LIBPATH} \
    -set prefix ${PREFIX}

echo "Formatting exonic region annotations"
sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/format_fusion_annotations.sas \
    -work $SASWORK \
    -set libpath ${LIBPATH} \
    -set prefix ${PREFIX}

echo "Format exon fragment annotations"
sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/format_exon_fragment_annotations.sas \
    -work $SASWORK \
    -set libpath ${LIBPATH} \
    -set prefix ${PREFIX}

echo "Flag exonic regions as constitutive, common, alternative"
sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/fusions_flag_constitutive_common_alt.sas \
    -work $SASWORK \
    -set libpath ${LIBPATH} \
    -set prefix ${PREFIX}

echo "Flag exon fragments as constitutive, common, alternative"
sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/fragments_flag_constitutive_common_unique.sas \
    -work $SASWORK \
    -set libpath ${LIBPATH} \
    -set prefix ${PREFIX}

echo "Make exon fragment and exonic region BED files for coverage"

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/fragments_export_bed_file.sas \
    -work $SASWORK \
    -set libpath ${LIBPATH} \
    -set prefix ${PREFIX} \
    -set outPath ${OUTPUT}

# Import logical junctions into SAS
echo "Importing all possible, logical junctions..."

### External variables to use in SAS program
INPUTFILE=${OUTPUT}/${PREFIX}_logical_junctions.bed
CHRSIZE=$(cut -f1 ${INPUTFILE} | wc -L)
EVENTIDSIZE=$(cut -f4 ${INPUTFILE} | wc -L)
BLOCKSIZE=$(cut -f11 ${INPUTFILE} | wc -L)
BLOCKSTARTS=$(cut -f12 ${INPUTFILE} | wc -L)
EXONSIZE=$(expr ${EVENTIDSIZE} / 2 + 1)

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/import_junctions_and_format.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set inputfile ${INPUTFILE} \
    -set chrsize ${CHRSIZE} \
    -set eventsize ${EVENTIDSIZE} \
    -set blocksize ${BLOCKSIZE} \
    -set blockstarts ${BLOCKSTARTS} \
    -set exonsize ${EXONSIZE}

# Import exon info
echo "Importing exon info for junction annotations..."

### External variables to use in SAS program
INPUTFILE=${OUTPUT}/${PREFIX}_exons.csv
CHRSIZE=$(cut -d, -f1 ${INPUTFILE} | wc -L)
EXONSIZE=$(cut -d, -f5 ${INPUTFILE} | wc -L)
XSCRIPTSIZE=$(cut -d, -f6 ${INPUTFILE} | wc -L)
GENESIZE=$(cut -d, -f7 ${INPUTFILE} | wc -L)
EVENTSIZE=$(expr ${EVENT_SIZE} \* 2)

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/import_exons_and_format.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set inputfile ${INPUTFILE} \
    -set chrsize ${CHRSIZE} \
    -set exonsize ${EXONSIZE} \
    -set xscriptsize ${XSCRIPTSIZE} \
    -set genesize ${GENESIZE} \
    -set eventsize ${EVENTSIZE}

# Import transcript junctions
echo "Importing junctions from transcripts..."

### External variables to use in SAS program
INPUTFILE=${OUTPUT}/${PREFIX}_transcript_junctions.csv
JUNCSIZE=$(cut -d, -f1 ${INPUTFILE} | wc -L)
JUNCCOORD=$(cut -d, -f2 ${INPUTFILE} | wc -L)
XSCRIPTSIZE=$(cut -d, -f3 ${INPUTFILE} | wc -L)
GENESIZE=$(cut -d, -f4 ${INPUTFILE} | wc -L)

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/import_junctions_from_transcripts.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set inputfile ${INPUTFILE} \
    -set juncsize ${JUNCSIZE} \
    -set junccoord ${JUNCCOORD} \
    -set xscriptsize ${XSCRIPTSIZE} \
    -set genesize ${GENESIZE}

# Import exon-skipping annotations
echo "Importing annotations for exon-skipping junctions..."

### External variables to use in SAS program
INPUTFILE1=${OUTPUT}/${PREFIX}_exon_skipping_junctions.csv
INPUTFILE2=${OUTPUT}/${PREFIX}_skipped_exons_list.csv
JUNCSIZE1=$(cut -d, -f1 ${INPUTFILE1} | wc -L)
SKIPPED1=$(cut -d, -f4 ${INPUTFILE1} | wc -L)
JUNCSIZE2=$(cut -d, -f1 ${INPUTFILE2} | wc -L)
SKIPPED2=$(cut -d, -f3 ${INPUTFILE2} | wc -L)

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/import_skipped_exon_annotations.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set inputfileA ${INPUTFILE1} \
    -set inputfileB ${INPUTFILE2} \
    -set juncsizeA ${JUNCSIZE1} \
    -set juncsizeB ${JUNCSIZE2} \
    -set skippedA ${SKIPPED1} \
    -set skippedB ${SKIPPED2}

# Import intron retention events
echo "Importing border junctions..."


### External variables to use in SAS program
INPUTFILE1=${OUTPUT}/${PREFIX}_intron_retention.csv
INPUTFILE2=${OUTPUT}/${PREFIX}_intron_retention.bed


GENESIZE=$(cut -d, -f1 ${INPUTFILE1} | wc -L)
EVENTIDSIZE=$(cut -d, -f2 ${INPUTFILE1} | wc -L)
CHRSIZE=$(cut -d, -f3 ${INPUTFILE1} | wc -L)
EXONSIZE=$(cut -d, -f6 ${INPUTFILE1} | wc -L)
EXONCATSIZE=$(cut -d, -f7 ${INPUTFILE1} | wc -L)
BLOCKSIZE=$(cut -f9 ${INPUTFILE2} | wc -L)
BLOCKSTARTS=$(cut -f10 ${INPUTFILE2} | wc -L)

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/import_intron_retention_events.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set inputfileA ${INPUTFILE1} \
    -set inputfileB ${INPUTFILE2} \
    -set genesize ${GENESIZE} \
    -set eventsize ${EVENTIDSIZE} \
    -set chrsize ${CHRSIZE} \
    -set splitsize ${EVENT_SIZE} \
    -set exonsize ${EXONSIZE} \
    -set exoncatsize ${EXONCATSIZE} \
    -set blocksize ${BLOCKSIZE} \
    -set blockstart ${BLOCKSTARTS}

# Flag junctions if derivable from transcripts
echo "Annotating logical junctions as transcript-derived junctions..."

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/flag_annotated_junctions.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH}

# Flag junctions if exon-skipping junctions
echo "Annotating logical junctions as exon-skipping junctions..."

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/flag_exonskip_junctions.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH}

# Concatenate exon-exon junctions and border junctions
echo "Stacking junctions..."

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/cat_junctions_and_ir.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set xscriptsize ${XSCRIPTSIZE}

# Add exon information to AS events
echo "Adding exon information to junctions..."

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/add_exon_info.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set genesize ${GENESIZE}

# Collapse duplicate events
echo "Collapsing duplicate junctions..."

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/collapse_duplicates.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set readSize ${READ_LENGTH} \
    -set eventsize ${EVENTIDSIZE} \
    -set genesize ${GENESIZE}

# Format and export splicing catalog
echo "Formatting annotation catalog and exporting..."

EVENTSIZE=$(expr ${EVENT_SIZE} \* 2)
OUTFILE=${OUTPUT}/${PREFIX}_splicing_catalog_${EVENTSIZE}bp.csv

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/format_splicing_annotations.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set outPath ${OUTFILE} 

###########################################################################################
#
# SECTION 3: CREATE REFERENCE SEQEUNCES AND BED FILES FOR COVERAGE
#
###########################################################################################

# Export a BED file for extracting sequences
echo "Exporting BED file of junctions..."

EVENTSIZE=$(expr ${EVENT_SIZE} \* 2)
OUTFILE=${OUTPUT}/${PREFIX}_splicing_catalog_${EVENTSIZE}bp.bed

sas -memsize 20G \
    -sysin ${SCRIPTS}/sas_programs/catalogue2BED.sas \
    -work ${SASWORK} \
    -set libpath ${LIBPATH} \
    -set outPath ${OUTFILE}

# Extract sequences
echo "Extracting reference FASTA file sequences..."

bedtools getfasta -fi ${FASTA} \
                  -bed ${OUTPUT}/${PREFIX}_splicing_catalog_${EVENTSIZE}bp.bed \
                  -fo ${OUTPUT}/${PREFIX}_splicing_catalog_${EVENTSIZE}bp.fa \
                  -name -s -split

# Make BED file for calculating coverage
echo "Making coverage BED file for splicing event FASTA sequences..."

python ${SCRIPTS}/fasta2bed.py --fasta ${OUTPUT}/${PREFIX}_splicing_catalog_${EVENTSIZE}bp.fa --out ${OUTPUT}/${PREFIX}_splicing_catalog_${EVENTSIZE}bp_for_coverage.bed

echo "Script done!"

