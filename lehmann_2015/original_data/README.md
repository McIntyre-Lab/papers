# Description of the CEGS F1-hybrid population

Raw RNA-seq data from the F1-hybrids have been deposited at NCBI as BioProject
PRJNA281652.

http://www.ncbi.nlm.nih.gov/bioproject/PRJNA281652/

We also rely heavily on:

PRJNA36679: Illumina WGS for Raleigh by DGRP
PRJNA74721: Illumina WGS for Winters

## Drosophila culture

The DGRP strains and the w[1118] strain were obtained from the Bloomington
Drosophila Stock Center. The Winter’s lines were collected and maintained in
the Nuzhdin laboratory. The w[1118] strain is the sequenced and isogenized
version (w[1118]-iso; 2-iso; 3-iso)[2, 3]. Culture and maintenance of
fly stocks were as in[1], with minor modifications.  The media
recipe was 150 grams sucrose, 150 grams yeast and 2 grams agar per liter of
water brought to full boil. Prior to dispensing media into vials, 15 ml of
tegosept stock solution (193 g/l in 95% ethanol) and 3 ml of propionic acid was
added per liter. The DGRP lines, Winter’s lines and the w[1118] strain were
cured of Wolbachia by three generations culture with doxycycline [4],
followed by more than 3 generations of culture in absence of doxycycline prior
to use.  Fly population densities were controlled when culturing the flies by
collecting eggs from population cages and pipetting equal number of eggs (32ul,
or ~300 eggs) into culture bottles. Virgins of the DGRP and Winter’s strains
were crossed to w[1118] males to generate hybrid progeny (cross #1). The
hybrid flies were collected on 15% yeast/sucrose media [5]. On day
5, virgin females were separated into two groups and one group was exposed to
w[1118] males (cross #2) for 24 hours. Both groups (mated and virgin) flies
were anesthetized using humidified CO2 gas, and in the mated group the females
were separated from the males, whereas the virgin controls were mock-sorted
using humidified CO2 gas. All flies were placed in plastic snap-cap tubes and
the tubes placed in liquid nitrogen. 

## DNA isolation

DNA was extracted from whole-body female flies using Qiagen’s DNeasy Blood and
Tissue Kit (Qiagen) and sheared to a fragment length of ~300 bp using the
Covaris S2 (Covaris). Subsequent library preparation was done according to
standard Illumina protocols.

## Mapping and SNP calling

Illumina sequences from two previous Drosophila genome resequencing projects
were downloaded from the NCBI Short Read Archive. These projects sequenced
genomic DNA from Drosophila sampled from Raleigh, NC (PRJNA36679) and Winters,
CA (PRJNA74721). For each line, the raw sequencing reads were trimmed by
quality using the SolexaQA package (ver. 1.12) with default parameters and all
trimmed reads less than 25 bp were discarded [14]. The quality trimmed reads
were then mapped to the *D. melanogaster* reference genome (FlyBase version 5.41)
using Bowtie 2 (ver. beta4) with the “very-sensitive” and “–N=1” parameters
[15]. Following mapping, the GATK (ver. 1.1-23)[16] IndelRealigner tool was
used to perform local realignments around indels and PCR and optical duplicates
were identified with the MarkDuplicates tool in the Picard package
(http://picard.sourceforge.net). SNP variants were identified in all lines
simultaneously using the GATK UnifiedGenotyper tool (ver. 2.1-8) with all
parameters set to the recommended default values[17]. In all, 216 Drosophila
lines were used for SNP calling. These lines included 179 lines from Raleigh,
NC (part of the DGRP project), 36 lines from Winters, CA, and the w[1118] line.

## SNP filtering

SNP calls were further filtered using VCFtools (version 1.1)[18] to remove all
sites with a genotyping rate less than 90%, sites with less than 2X coverage,
and all sites with more than two alleles. All of the lines included in this
study have been subjected to many generations of inbreeding prior to
sequencing. While there is always a possibility of some polymorphic positions
due to residual heterozygosity and/or new mutations, most loci are expected to
be homozygous. SNPs identified as heterozygous were filtered out and not
included in subsequent analysis. Variants with fewer than 10 lines having
evidence for a particular polymorphism were also excluded from analysis. The
2.1 million SNP calls for the Raleigh lines were compared to release 2.0 [3].
There were only minor differences in the SNP calls (95% were identical) and
differences were largely due to rare variants, some of which were present in
release 2.0 and absent in our calls and vice versa. This low level of
discrepancy is not surprising given the contributions of the additional Winters
lines to variant calling. 

## mRNA 

For each genotype and condition (mated or virgin), the flies were pooled to
create three replicates of approximately 50 flies each. Heads were separated
from bodies by shaking frozen flies, and the heads were sorted from the bodies
using the dissecting microscope and placed immediately in trizol. Isolation of
mRNA from the total RNA pool was done using the Dynabeads® mRNA purification
kit (Invitrogen Dynal AS, Oslo, Norway). Final elution of pure mRNA was done in
15 µl of Tris-HCl (pH 7.5). 6µl were kept for long-term storage, and the
remaining 9 µl were used for subsequent library construction. Libraries were
made in a 96-well plate format following [6]. Indexing of the samples was done
using a combination of eight Illumina indexes (ATCACG, CGATGT, TTAGGC, TGACCA,
ACAGTG, GCCAAT, CAGATC, ACTTGA) and twelve 4-nucleotides barcodes (ACTG, ATGC,
AGCT, TACG, TCGA, TGAC, CAGT, CTAG, CGTA, GATC, GCAT, GTCA) ligated at both
ends of the fragment. In this way, each of the 96 samples in the plate has a
unique combination of index-barcode. All samples in a plate were pooled
together in a single tube for sequencing. These pooled samples were visualized
by Agilent Bioanalyzer and quantified using the KapaBiosystems Library
Quantification Kit, according to manufacturer’s instructions. The pooled
samples were loaded into an Illumina flow cell v.3 at a concentration of 16 pM
and run on the HiSeq 2000 for 2x100 cycles. Library quality control and
sequencing was performed at the USC NCCC Epigenome Center’s Data Production
Facility.  Certain samples were sequenced as Single End and others as Paired
End, reflecting the change in cost over the term of this project. Each sample
was sequenced over multiple lanes [44].

Illumina sequences were demultiplexed into each of the individual indices using
CASAVA and allowing for one mismatch. As additional quality control the
internal barcodes that were present on both ends of the fragment were verified.
In rare cases, one end did not match to any of the existing barcodes, but the
other internal barcode was a perfect match and the read was assigned to the
sample with the barcode that was a perfect match. All Illumina data are
available on the NCBI Short Read Archive under the BioProject PRJNA74721.

Distinct reads (those with no exact sequence duplicates) were first aligned to
a set of all logical junctions created from FlyBase annotations using Bowtie
(v0.12.9, -v3 -m1)[7].  Distinct reads that did not uniquely align to a
junction were aligned to a non-redundant set of exonic regions. The 63,706
exonic regions of Flybase version r5.51 were compared, and 762 of these regions
had an exact sequence duplication with at least one other region. These
duplicated regions were reduced to a single representative copy, with the
multiple gene annotation noted for a total of 63,181 distinct exonic regions.
Note that this will not eliminate map bias in areas of sequence similarity or
genome ambiguity [8]. Reads were aligned uniquely to this non-redundant exon
reference using Bowtie (v0.12.9, with options -v3 -m1). Unaligned reads were
then aligned using LAST which allows for soft trimming and thus minimizes bias
(v193, -l25) [10]. Initial Coverage was assessed at the exonic region level by
calculating the average number of reads per nucleotide coverage (APN) or the
number of reads divided by the length of the region [11]. The performance of
this sequential alignment process was compared to alignments generated using
BWA-mem.  In this comparison the whole genome was used a reference. No
significant differences were observed between the two approaches in terms of
number of reads in exonic regions.  

Samples with fewer than 50% of all exonic regions with at least one read were
not considered further (n=104). Next exonic regions were filtered based on
coverage. An exonic region was considered expressed if >50% of all replicates
had at least one read for each genotype by mating status. A total of 35,676
exonic regions (virgin) and 37,442 exonic regions (mated) were retained with
this criteria. UQ normalization [13] was selected. Normalization factors were
calculated separately for mated and virgin environments.  After normalization
the natural log was taken and the genotype*mating status was centered using the
upper quartile [13].   
