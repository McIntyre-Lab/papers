#! /bin/bash
# Marty McCrory, 2011-12-29

work=/scratch/hpc/mccrory/fru_network

for gene in dsx "fl(2)d" fru her ix msl-2 snf Sxl tra2 tra
do
	for sample in AH_BerF AH_BerM AH_CS AH_dsxD AH_dsxNullF AH_dsxNullM AH_FruP14-440 AH_FruW12-ChaM5 AH_Female_FruM-A AH_Male_FruM-A AH_Female_FruM-B AH_Male_FruM-B AH_Female_FruM-C AH_Male_FruM-C AH_CSFemale
	do
		date1=0; barcode1=0; lane1=0
		date2=0; barcode2=0; lane2=0
		date3=0; barcode3=0; lane3=0
		date4=0; barcode4=0; lane4=0
		
		case $sample in
			"AH_BerF")
				date1=2011-05-03; barcode1=TGACCA; lane1=5
				date2=2011-05-03; barcode2=ACAGTG; lane2=3
				date3=2011-05-03; barcode3=GCCAAT; lane3=2
				;;

			"AH_BerM")
				date1=2011-05-03; barcode1=ATCACG; lane1=1
				date2=2011-05-03; barcode2=CGATGT; lane2=6
				date3=2011-05-03; barcode3=TTAGGC; lane3=5
				;;

			"AH_CS")
				date1=2011-05-03; barcode1=ACTTGA; lane1=1
				date2=2011-05-03; barcode2=ACTTGA; lane2=2
				date3=2011-05-03; barcode3=ACTTGA; lane3=3
				date4=2011-05-03; barcode4=ACTTGA; lane4=5
				;;

			"AH_dsxD")
				date1=2011-05-03; barcode1=CAGATC; lane1=1
				date2=2011-05-03; barcode2=ACTTGA; lane2=6
				date3=2011-05-03; barcode3=GATCAG; lane3=5
				;;

			"AH_dsxNullF")
				date1=2011-05-03; barcode1=TAGCTT; lane1=6
				date2=2011-05-03; barcode2=ATCACG; lane2=2
				date3=2011-05-03; barcode3=CGATGT; lane3=5
				;;

			"AH_dsxNullM")
				date1=2011-05-03; barcode1=TTAGGC; lane1=6
				date2=2011-05-03; barcode2=TGACCA; lane2=3
				date3=2011-05-03; barcode3=ACAGTG; lane3=6
				;;

			"AH_FruP14-440")
				date1=2011-05-03; barcode1=GATCAG; lane1=1
				date2=2011-05-03; barcode2=GATCAG; lane2=2
				date3=2011-05-03; barcode3=GATCAG; lane3=3
				;;

			"AH_FruW12-ChaM5")
				date1=2011-05-03; barcode1=TAGCTT; lane1=1
				date2=2011-05-03; barcode2=TAGCTT; lane2=2
				date3=2011-05-03; barcode3=TAGCTT; lane3=3
				date4=2011-05-03; barcode4=TAGCTT; lane4=7
				;;

			"AH_Female_FruM-A")
				date1=2011-07-05; barcode1=TGACCA; lane1=1
				date2=2011-07-05; barcode2=ACAGTG; lane2=2
				date3=2011-07-05; barcode3=GCCAAT; lane3=3
				;;

			"AH_Male_FruM-A")
				date1=2011-07-05; barcode1=ATCACG; lane1=1
				date2=2011-07-05; barcode2=CGATGT; lane2=2
				date3=2011-07-05; barcode3=TTAGGC; lane3=3
				;;

			"AH_Female_FruM-B")
				date1=2011-07-05; barcode1=TGACCA; lane1=2
				date2=2011-07-05; barcode2=ACAGTG; lane2=3
				date3=2011-07-05; barcode3=GCCAAT; lane3=1
				;;

			"AH_Male_FruM-B")
				date1=2011-07-05; barcode1=ATCACG; lane1=2
				date2=2011-07-05; barcode2=CGATGT; lane2=3
				date3=2011-07-05; barcode3=TTAGGC; lane3=1
				;;

			"AH_Female_FruM-C")
				date1=2011-07-05; barcode1=TAGCTT; lane1=5
				date2=2011-07-05; barcode2=ATCACG; lane2=6
				date3=2011-07-05; barcode3=CGATGT; lane3=5
				;;

			"AH_Male_FruM-C")
				date1=2011-07-05; barcode1=CAGATC; lane1=1
				date2=2011-07-05; barcode2=ACTTGA; lane2=2
				date3=2011-07-05; barcode3=GATCAG; lane3=3
				;;

			"AH_CSFemale")
				date1=2011-07-05; barcode1=CAGATC; lane1=7
				date2=2011-07-05; barcode2=ACTTGA; lane2=7
				date3=2011-07-05; barcode3=GATCAG; lane3=8
				;;	
		esac
	
	
		sample1=${work}/gene_pileups/cut_cols/${date1}_${lane1}_${barcode1}_${gene}.tsv
		sample2=${work}/gene_pileups/cut_cols/${date2}_${lane2}_${barcode2}_${gene}.tsv
		sample3=${work}/gene_pileups/cut_cols/${date3}_${lane3}_${barcode3}_${gene}.tsv
		sample4=${work}/gene_pileups/cut_cols/${date4}_${lane4}_${barcode4}_${gene}.tsv
		
		output=${work}/gene_pileups/cut_cols/adult_heads_combined_samples/${sample}_${gene}.tsv
	
	
		cat $sample1 $sample2 $sample3 $sample4 2>/dev/null | perl ${work}/scripts/combine_cut_cols_pileup_samples.pl >$output
	done
done

