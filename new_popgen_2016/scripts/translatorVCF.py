#!/usr/bin/env python

#Standart Libraries
import re
import os


vcf_file = open('no_missing_biall_chr_nolab_nomel.recode.vcf','r')
outdir = ''

for line in vcf_file:
	#skiping comentaries but not header
	if line.startswith('##'):continue

	#Taking out newline and spliting it for tabs into a list
	line = line.strip('\n').split('\t')

	#Catching chromosome for filename
	filename = line[0]

	#Removing 0/0:4,0:4:12:0,12,180 everything afther the first ':' 
	line = [element.split(':')[0] if re.match('[0-9]/.+',element) else element  for element in line]

	#Translating 0/0 to reference, 0/1 to N and 1/1 to alt.
	line = [line[3] if element == '0/0' else line[4] if element == '1/1' else 'N' if element == '0/1' else element for element in line]

	#Getting the new header
	if line[0].startswith('#'):
		header = line[1:2]+line[9:]
		header[0] = '#'+header[0]
		header = ','.join(header)
		line=''

	#Opening output file per chromosome
	else: 
		output = open(filename+'.h12','a')
		line = line[1:2]+line[9:]
		line = ','.join(line)

	#Writting out 
	if line and header:
		print >> output,header
		print >> output,line
		header = ''
	elif line:
		print >> output,line
	else:
		pass
	#closing the file if open
	if line:
		output.close