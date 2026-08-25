## Orthologous Transcripts and Alternative Splicing in Drosophila
This paper (bankole_sex_specific_splicing_2026) contains the documentation and scripts 
used to identify orthologous transcripts and characterize alternative splicing in Drosophila species

The analysis uses long read RNAseq data from male and female head tissue to identify transcript models, 
quantify transcripts and characterize alternative splicing across multiple Drosophila species.

## Species:
Drosophila melanogaster (dmel6)
Drosophila simulans (dsim2)
Drosophila sechellia (dser1)
Drosophila yakuba (dyak2)
Drosophila santomea (dsan1)

## Overview:
Transcript homology is inferred using reciprocal liftover of
transcript models between species. Orthologous transcript groups are
then identified using network-based analyses of the reciprocal
mapping relationships.

The resulting transcript models and orthology assignments are
validated using long read RNA-seq data from male and female head tissue.


## Primary analysis steps:
1. Create fiveSpecies transcript model annotations based on Unique Junction Chains (UJCs)
2. Generate UJCs and Exon Region Patterns (ERPs) from long read RNAseq data
3. Perform SQANTI3 QC and SQANTIreads QC analyses
4. Prepare data for UJC and ERP analysis
5. Normalize gene level expression
6. Generate gene level summaries
7. Identify and characterize alternative splicing
8. Generate component-level summaries


## Organization 
bash_scripts == Bash / submission scripts to run analyses  
design_files == Sample and analysis design files
documentation == Supporting documentation
sas_programs == SAS programs for statistical analysis
scripts == Python/R analysis scripts
sra_info == SRA metadata and related information

### Unique Junction Combination (UJC)
A UJC represents a unique combination of splice junctions within a
transcript model. UJCs are used as a primary unit for transcript
quantification and cross-species transcript comparisons. We use a 64-character
hexadecimal hash (jxnHash) for computational efficiency and data 
management.

### Transcript Homology
Transcript homology is inferred through reciprocal liftover between
species. Transcript models that map reciprocally between species are
used to construct networks of corresponding transcript models.

### Orthology
Orthologous transcript groups are identified from the network of
reciprocal transcript mappings. Each connected component represents a
group of transcript models with inferred evolutionary correspondence.

### Alternative Splicing
Alternative splicing is characterized using UJC and ERP-level
expression data to identify differences in transcript structure and
splicing between species and between males and females.

## Reproducibility
The repository contains the scripts, design files, documentation, and
statistical programs used in the analysis. Individual analysis steps
are documented in the documentation directory.  

## Additional data availability
The fiveSpecies full annotation files, supporting files, supplemental files, tables, figures
are available on zenodo:   [link!!!]

