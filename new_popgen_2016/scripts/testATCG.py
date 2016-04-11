#!/usr/bin/env python
#Standart Libraries
import re

vcf_file = open('no_missing_biall_chr_nolab_nomel.recode.vcf','r')
noATCGfound=0
for line in vcf_file:
	if line.startswith('#'):continue
	line = line.strip('\n')
	line = line.split('\t')
	line = line[3]+line[4]
	line = line.replace('A','')
	line = line.replace('T','')
	line = line.replace('C','')
	line = line.replace('G','')
	noATCGfound=noATCGfound+len(line)

print "No ATCG letters found: \t",noATCGfound