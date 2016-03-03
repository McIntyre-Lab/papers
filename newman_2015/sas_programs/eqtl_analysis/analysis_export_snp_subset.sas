/* Export subsets of SNPs for calculating LD */

/* Set libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Extract list of credible SNPs from the list of genotyped SNPs */
/* Need to make UCSC-coordinate style IDs for the credible SNPs, as some ImmunoChip SNPs don't have rs# */

data credible_snps;
  set eqtl.supptable1_snp_list;
  length ucsc_snp_id $15.;
  ucsc_snp_id=cats(chr,':',position);
  keep snp_id ucsc_snp_id;
run;

data genotyped_snps;
   set eqtl.genotyped_snps;
   keep snp_id;
run;

/* Extract credible SNPs from all genotyped SNP list - rs# */

proc sort data=credible_snps nodup;
   by snp_id;
proc sort data=genotyped_snps nodup;
   by snp_id;
run;

data credible_snp_list no_cred_rs no_genotypes;
   merge genotyped_snps (in=in1) credible_snps (in=in2);
   by snp_id;
   if in1 and in2 then output credible_snp_list;
   else if in1 then output no_cred_rs;
   else output no_genotypes;
run;

* 2027 possible SNPs from Onengut-Gumuscu, 2015 paper;
* 164643 SNPs with genotype information;
* 1247 O-G SNPs genotyped AND have an rs#;
* 780 O-G SNPs that either have no genotypes, or no rs#;
* 163396 genotyped SNPs remaining;

/* Extract credible SNPs from remaining genotyped SNP list - UCSC ID */

data cred_snps_ucsc;
   set no_genotypes;
   rename snp_id=rs_snp_id ucsc_snp_id=snp_id;
run;

proc sort data=cred_snps_ucsc;
   by snp_id;
proc sort data=no_cred_rs;
   by snp_id;
run;

data credible_snp_list_ucsc no_cred_snp no_genotype_info;
   merge no_cred_rs (in=in1) cred_snps_ucsc (in=in2);
   by snp_id;
   if in1 and in2 then output credible_snp_list_ucsc;
   else if in1 then output no_cred_snp;
   else output no_genotype_info;
run;

* 780 possible UCSC-style SNPs from Onengut-Gumuscu, 2015 paper;
* 163396 remaining SNPs with genotype information;
* 7 UCSC-style O-G SNPs genotyped AND have a UCSC ID;
* 773 O-G SNPs that have no genotypes;
* 163389 genotyped SNPs remaining;

/* Quick check: rs## might be different! */

data check_snps_wo_gts;
   set no_genotype_info;
   keep snp_id rs_snp_id;
run;

data genotyped_snp_ucsc_id;
   set eqtl.genotyped_snps;
   length snp_coord_id $15.;
   snp_coord_id=cats('chr',chr,':',pos);
   keep snp_id snp_coord_id;
   rename snp_id=snp_id_rs snp_coord_id=snp_id;
run;

proc sort data=check_snps_wo_gts;
   by snp_id;
proc sort data=genotyped_snp_ucsc_id;
   by snp_id;
run;

data genotyped_snp_bad_id genotyped_snp_good_id no_genos;
   merge genotyped_snp_ucsc_id (in=in1) check_snps_wo_gts (in=in2);
   by snp_id;
   if in1 and in2 then output genotyped_snp_bad_id;
   else if in1 then output genotyped_snp_good_id;
   else output no_genos;
run;

*1 O-G SNP with a different ID than on ImmunoChip;
*772 O-G SNPs without any genotype info;

/* I am going to make an "ImmunoChip ID"-to-"Credible SNP ID" index dataset. This is so that there are no issues merging between different IDs */

/* I will have:
1. ImmunoChip ID
2. O-G ID
3. UCSC Coordinate ID
4. "universal ID", which will be based on the ImmunoChip ID
*/

data credible_snp_list_2;
   set credible_snp_list;
   length ic_snp_id $15.;
   length onengut_snp_id $15.;
   onengut_snp_id=snp_id;
   ic_snp_id=snp_id;
run;

data credible_snp_list_ucsc_2;
   set credible_snp_list_ucsc;
   length ic_snp_id $15.;
   length onengut_snp_id $15.;
   length ucsc_snp_id $15.;
   onengut_snp_id=rs_snp_id;
   ic_snp_id=snp_id;
   ucsc_snp_id=snp_id;
   keep snp_id ucsc_snp_id ic_snp_id onengut_snp_id;
run;

data genotyped_snp_bad_id_2;
   set genotyped_snp_bad_id;
   length uni_snp_id $15.;
   length ic_snp_id $15.;
   length onengut_snp_id $15.;
   length ucsc_snp_id $15.;
   uni_snp_id=snp_id_rs;
   ic_snp_id=snp_id_rs;
   ucsc_snp_id=snp_id;
   onengut_snp_id=rs_snp_id;
   keep uni_snp_id ucsc_snp_id ic_snp_id onengut_snp_id;
   rename uni_snp_id=snp_id;
run;

data onengut_snps_no_geno;
   set no_genos;
   length uni_snp_id $15.;
   length ic_snp_id $15.;
   length onengut_snp_id $15.;
   length ucsc_snp_id $15.;

   uni_snp_id='';
   ic_snp_id='';
   ucsc_snp_id=snp_id;
   onengut_snp_id=rs_snp_id;
   keep uni_snp_id ucsc_snp_id ic_snp_id onengut_snp_id;
   rename uni_snp_id=snp_id;
run;

data eqtl.ic_to_onengut_snp_index;
   set credible_snp_list_2 credible_snp_list_ucsc_2 genotyped_snp_bad_id_2 onengut_snps_no_geno;
run;

/* Extract eQTL SNPs from all genotyped SNP list */


data eqtl_tested_snps;
   set eqtl.eqtl_results_summary_table;
   keep snp_id;
run;

* remove duplicated SNP ids;
proc sort data=eqtl_tested_snps nodup;
   by snp_id;
run;

proc sort data=genotyped_snps;
   by snp_id;
run;

data eqtl_tested_snp_list no_test no_geno_oops;
   merge genotyped_snps (in=in1) eqtl_tested_snps (in=in2);
   by snp_id;
   if in1 and in2 then output eqtl_tested_snp_list;
   else if in1 then output no_test;
   else output no_geno_oops;
run;

*164643 SNPs with genotypes, 11368 SNPs tested for eQTLs;
*11368 eQTL-SNPs in genotyped SNP list, 153275 genotyped SNPs not tested for eQTLs;
*0 eQTL-SNPs without genotypes;

/* Stack lists and drop duplicates */

data ic_assoc_snps;
   set eqtl.ic_to_onengut_snp_index;
   where snp_id ne '';
   keep snp_id;
run;

data snps_to_subset;
   set ic_assoc_snps eqtl_tested_snp_list;
run;

proc sort data=snps_to_subset nodup;
   by snp_id;
run;

/* Need to split the eQTL tested SNP list into 2 lists for PLINK to not give weird results
   Going to arbitrarily group the SNPs into two chromosome groups: chr1-chr9 and chr10-chrY
   as this seems to fix the issue */

data eqtl_snp2chrom;
  set eqtl.snp_data_w_info;
  keep snp_id chr;
run;

proc sort data=eqtl_snp2chrom nodup;
   by snp_id;
proc sort data=eqtl_tested_snp_list;
   by snp_id;
run;

data eqtl_tested_snps_w_chr;
    merge eqtl_snp2chrom (in=in1) eqtl_tested_snp_list (in=in2);
   by snp_id;
   if in1 and in2;
run;

proc sort data=eqtl_tested_snps_w_chr;
   by chr;
run;


data eqtl_tested_snp_list_1 eqtl_tested_snp_list_2;
   set eqtl_tested_snps_w_chr;
   if chr*1 lt 10 then output eqtl_tested_snp_list_1;
   else output eqtl_tested_snp_list_2;
   keep snp_id;
run;

/* Export SNP lists */

proc export data=snps_to_subset outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/subset_snp_list.txt'
   dbms=tab replace; putnames=no; run;

proc export data=eqtl_tested_snp_list_1 outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_snp_list_1.txt'
   dbms=tab replace; putnames=no; run;

proc export data=eqtl_tested_snp_list_2 outfile='/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eqtl_snp_list_2.txt'
   dbms=tab replace; putnames=no; run;

