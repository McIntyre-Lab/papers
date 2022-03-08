#!/bin/sh
#SBATCH --mail-user=adalena.nanni@ufl.edu
#SBATCH --job-name=bwa
#SBATCH --account=mcintyre
#SBATCH --qos=mcintyre-b
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue
#SBATCH -o /blue/mcintyre/share/maize_ainsworth/scripts/rnaseq/SLURM_LOGS/bwa.%A.out
#SBATCH -t 04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6gb

### create BWA reference from Mo17 YAN and Mo17 CAU genomes

#Load modules
module load bwa/0.7.15

#Set directories
PROJ=/blue/mcintyre/share/maize_ainsworth

##### Mo17 YAN Genome #####

REFGYAN=/blue/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN

#### create bwa ref for mo17 yan genome
#BWA-mem references  -p is the prefix for the output database,
#    -a is the indexing algorithm ('bwtsw' is for ref>2G, 'is' for ref<2G).
bwa index -p $REFGYAN/Zm-Mo17-REFERENCE-YAN-1.0_BWA \
    -a  bwtsw $REFGYAN/Zm-Mo17-REFERENCE-YAN-1.0.fasta


##### Mo17 CAU Genome #####

REFGCAU=/blue/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_CAU

#### create bwa ref for mo17 cau genome
#BWA-mem references  -p is the prefix for the output database,
#    -a is the indexing algorithm ('bwtsw' is for ref>2G, 'is' for ref<2G).
bwa index -p $REFGCAU/Zm-Mo17-REFERENCE-CAU-1.0_BWA \
    -a  bwtsw $REFGCAU/Zm-Mo17-REFERENCE-CAU-1.0.fsa
