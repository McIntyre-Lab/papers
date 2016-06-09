#!/bin/bash

tmpdir=`mktemp -d`

vcfdir=$HOME/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files
bed=$MCLAB/useful_dmel_data/flybase551/si_fusions/fb551_si_fusions.bed
vcfScript=$HOME/devel/mcscript/vcfCountSnps.py

for line in $(cat $MCLAB/cegs_ase_paper/design_files/CEGS_list_68_lines.txt); do
    $vcfScript --vcf $vcfdir/${line}_w11182${line}_UPD.vcf.gz --bed $bed --out $tmpdir/line.csv
    $vcfScript --vcf $vcfdir/w1118_w11182${line}_UPD.vcf.gz --bed $bed --out $tmpdir/tester.csv
    awk -v l=$line '{if (NR == 1){
                            gsub("regionID", "fusion_id", $0)
                            print "line," $0
                        }else{
                            print l "," $0
                        }
                    }' $tmpdir/line.csv > $tmpdir/${line}.csv

    awk -v l=w1118_$line '{if (NR == 1){
                                gsub("regionID", "fusion_id", $0)
                                print "line," $0
                            }else{
                                print l "," $0
                            }
                    }' $tmpdir/tester.csv > $tmpdir/w1118_${line}.csv
done
rm $tmpdir/line.csv
rm $tmpdir/tester.csv

/home/jfear/devel/mcscript/catTable.py -f $tmpdir/*.csv --header --odir /home/jfear/mclab/cegs_ase_paper/pipeline_output/ase_summary --oname snp_indel_cnt_by_fusion.csv
