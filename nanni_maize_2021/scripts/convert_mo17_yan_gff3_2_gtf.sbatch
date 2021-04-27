#!/bin/bash
#SBATCH --mail-user=ammorse@ufl.edu
#SBATCH --mail-type=FAIL,END
#SBATCH --job-name=gff2gtf
#SBATCH --output=/ufrc/mcintyre/share/maize_ainsworth/scripts/pacbio/SLURM_LOGS/gff2gtf_%A.out
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16gb
#SBATCH --cpus-per-task=1


### convert mo17 yan gff3 to gtf using gffread from cufflinks

ml cufflinks/2.2.1.1

INPUT=/ufrc/mcintyre/share/references/maize_mo17/mo17_genome_maizegdb_YAN

gffread -E -O -T $INPUT/Zm-Mo17-REFERENCE-YAN-1.0.mw2.gff3 -o $INPUT/Zm-Mo17-REFERENCE-YAN-1.0.mw2.gtf


