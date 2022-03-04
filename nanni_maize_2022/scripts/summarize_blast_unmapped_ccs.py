#!/usr/bin/env python

import pandas as pd
import os

def list_files(directory, extension, allow_ext=False):
    if allow_ext == True:
        return (f for f in os.listdir(directory) if f.endswith('.' + extension) or "."+extension+"." in f)
    else:
        return (f for f in os.listdir(directory) if f.endswith('.' + extension))

# Loop over *.tsv files in input directory to get all BLAST results
inDir = "~/mclab/SHARE/McIntyre_Lab/maize_ozone_FINAL/2018/PacBio/blast_unmapped_ccs"
blastDF = pd.DataFrame()
for file in list_files(inDir, "tsv"):
    tempDF = pd.read_csv(inDir + "/" + file, sep="\t")
    tempDF["sampleID"] = "_".join(file.split("_")[:-5])
    blastDF = pd.concat([blastDF, tempDF], ignore_index=True)

print(blastDF.groupby("sampleID")["qseqid"].nunique())
#    sampleID
#    113_c123_amb    352
#    120_c123_oz     354
#    19_mo17_amb     144
#    21-2_mo17_oz    233
#    21_mo17_oz      157
#    42_b73_amb      432
#    46_b73_oz       356
#    67_hp301_amb    456
#    70_hp301_oz     477
#    89_nc338_amb    466
#    96_nc338_oz     597
# Total number of reads with a blast hit (7637 total unmapped reads went into blast)
print(blastDF["qseqid"].nunique())
#    4024
# Get best hit for each read
bestHitDF = blastDF.groupby("qseqid")[
        [c for c in blastDF.columns if c != "qseqid"]
    ].first().reset_index()

# Get proportion of blast length to best hit (length/slength) of best hits
bestHitDF["length_prop_slen"] = bestHitDF["length"] / bestHitDF["slen"]

# Get distribution of length, evalue, mismatches, gaps, length_prop_slen across best hits
bestHitDF[["length", "evalue", "mismatch", "gaps", "length_prop_slen"]].describe().round(5)
#           length      evalue    mismatch        gaps  length_prop_slen
#count  4024.00000  4024.00000  4024.00000  4024.00000        4024.00000
#mean   1091.59990     0.00140    14.36854    39.44036           0.54252
#std     931.20902     0.01464    25.69250    58.56905           0.37210
#min      28.00000     0.00000     0.00000     0.00000           0.00000
#25%     262.00000     0.00000     2.00000     9.00000           0.21032
#50%     999.00000     0.00000     6.00000    19.00000           0.49306
#75%    1699.25000     0.00000    15.00000    42.00000           0.95615
#max    6792.00000     0.34000   311.00000   766.00000           1.21988

# Count how many best hits have 0 mismatches and 0 gaps
len(bestHitDF[bestHitDF["mismatch"] + bestHitDF["gaps"] == 0])
#   95
# Count how many best hits have 0 mismatches and 0 gaps with >0.5 length_prop_slen
len(bestHitDF[
        (bestHitDF["mismatch"] + bestHitDF["gaps"] == 0)&
        (bestHitDF["length_prop_slen"] > 0.5)])
#   2
print(bestHitDF[
        (bestHitDF["mismatch"] + bestHitDF["gaps"] == 0)&
        (bestHitDF["length_prop_slen"] > 0.5)
    ][["stitle","length"]].to_string(index=False))
#                                           stitle  length
#              Zea mays clone 400604 mRNA sequence     636
# Zea mays thioredoxin H-type (LOC100281021), mRNA     596

