#!/bin/bash 

## extract transcriptID, geneID and associated geneNames from GTF

INPUT=/home/ammorse/TB14/maize_gtf_files
OUTPUT=/home/ammorse/mclab/SHARE/McIntyre_Lab/useful_maize_info/RefGenV4

##GTF_extract --fields=transcript_id,gene_id,gene_name $INPUT/Zea_mays.B73_RefGen_v4.41.gtf | sort -u > $OUTPUT/B73_v4_transcriptIDs_geneIDs_geneNames_02amm.txt

GTF_extract --fields=gene_id,gene_name $INPUT/Zea_mays.B73_RefGen_v4.41.gtf | sort -u > $OUTPUT/B73_v4_geneIDs_geneNames_02amm.txt

