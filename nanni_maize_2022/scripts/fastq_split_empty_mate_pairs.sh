#!/bin/sh

## Takes in two arguments: R1.fq R2.fq
## Splits each file into *.no_empty_pairs.fq and *.with_empty_pairs.fq
## Output files are prefixed with input file names and places in same directory

## Get file prefixes
PREFIX1=$(dirname $1)/$(basename $(basename $1 .fq) .fastq)
PREFIX2=$(dirname $2)/$(basename $(basename $2 .fq) .fastq)


## Paste together R1 and R2 reads

## Check read lengths, if one is 0 then put in the 'with_empty_pair' files
paste $1 $2 | paste - - - - | \
    awk -F "\t" -v OFS="\n" '$3=="" || $4==""{print($1,$3,$5,$7)}' \
    > ${PREFIX1}.with_empty_pairs.fq
paste $1 $2 | paste - - - - | \
    awk -F "\t" -v OFS="\n" '$3=="" || $4==""{print($2,$4,$6,$8)}' \
    > ${PREFIX2}.with_empty_pairs.fq

## Check read lengths, if all reads present then put in the 'no_short_pair' files
paste $1 $2 | paste - - - - | \
    awk -F "\t" -v OFS="\n" '$3!="" && $4!=""{print($1,$3,$5,$7)}' \
    > ${PREFIX1}.no_empty_pairs.fq
paste $1 $2 | paste - - - - | \
    awk -F "\t" -v OFS="\n" '$3!="" && $4!=""{print($2,$4,$6,$8)}' \
    > ${PREFIX2}.no_empty_pairs.fq

## Check that line numbers match
CK1=$(cat ${PREFIX1}.with_empty_pairs.fq ${PREFIX1}.no_empty_pairs.fq | wc -l)
CK2=$(cat ${PREFIX2}.with_empty_pairs.fq ${PREFIX2}.no_empty_pairs.fq | wc -l)
if [[ $(cat $1 | wc -l) != ${CK1} ]]; then
    echo "
WARNING : Incorrect number of reads in split empty mate pairs of $1
"
else
    echo "Correctly split empty mate pairs of $1
"
fi
if [[ $(cat $2 | wc -l) != ${CK2} ]]; then
    echo "WARNING : Incorrect number of reads in split empty  mate pairs of $2"
else
    echo "Correctly split empty mate pairs of $2
"
fi
