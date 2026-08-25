## Orthologous Transcripts and Alternative Splicing in Drosophila
This paper directiory (bankole_2026) contains the documentation and scripts 
used to identify orthologous transcripts and characterize transcript usage in 
Drosophila species

The analysis uses long read RNAseq data from male and female head tissue to identify transcript models, 
quantify transcript expression and characterize transcript usage across multiple Drosophila species.

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
validated using long read RNAseq data from male and female head tissue.


## Primary analysis steps:
1. Create fiveSpecies transcript model annotations based on Unique Junction Chains (UJCs)
2. Generate UJCs and Exon Region Patterns (ERPs) from long read RNAseq data
3. Perform SQANTI3 QC and SQANTIreads QC analyses
4. Prepare data for UJC and ERP analysis
5. Normalize gene level expression
6. Generate transcript and gene level expression summaries
7. Identify and characterize alternative transcript usage


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
species. Transcript models that map between species are
used to construct networks of corresponding transcript models.

### Orthology
Orthologous transcript groups are identified from the network of
transcript mappings. Only components with exactly one transcript per 
species, either annotated or observed.

## Reproducibility
The repository contains the scripts, design files, documentation, and
statistical programs used in the analysis. Individual analysis steps
are documented in the documentation directory.  

## Additional data availability
The fiveSpecies full annotation files, supporting files, supplemental files, tables, figures
are available on zenodo:   [link!!!]
Data is available on SRA:  ######s

## IGV browser for visualizing transcript models
A separate repository (IGV Browser for papers/bankole_2026) containing the files used to 
visualize transcript models and transcript usage using Integrative Genomics Viewer (IGV).

The IGV tracks provide a visual representation of transcript models,
splice junctions, and long read RNAseq alignments.

- visualize transcript models across Drosophila species
- visualize transcript models identified from long read RNAseq
- compare transcript structures between males and females
- examine splice junction usage
- evaluate evidence supporting alternative splicing




