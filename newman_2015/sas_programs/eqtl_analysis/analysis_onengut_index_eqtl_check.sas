/* Set libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Check if any Onengut SNP group (index snp/region) is excluded from eQTL analysis */

* get the list of index snps for O-G SNPs that are in LD with an eQTL-tested SNP;
* For this we will assume r2>0.8 for high LD, or less than that for low LD;
data onengut_index_snps_high_ld onengut_index_snps_low_ld;
   set eqtl.onengut_snps_ld_eqtl_snps;
   if r2 ge 0.8 then output onengut_index_snps_high_ld;
   else output onengut_index_snps_low_ld;
   keep onengut_index_snp;
run;

* get the list of index snps for O-G SNPs that have no LD results;

data onengut_index_snps_no_ld;
   set eqtl.onengut_snps_ld_no_eqtl_snps;
   keep onengut_index_snp;
run;

* remove duplicates;

proc sort data=onengut_index_snps_high_ld nodup;
   by onengut_index_snp;
proc sort data=onengut_index_snps_low_ld nodup;
   by onengut_index_snp;
proc sort data=onengut_index_snps_no_ld nodup;
   by onengut_index_snp;
run;

data onengut_index_snp_eqtl_check;
  length index_group $4.;
  merge onengut_index_snps_high_ld (in=in1) onengut_index_snps_low_ld (in=in2) onengut_index_snps_no_ld (in=in3);
  by onengut_index_snp;
  if in1 then index_group='HIGH';
  else if in2 then index_group='LOW';
  else index_group='NONE';
run;

proc freq data=onengut_index_snp_eqtl_check;
   tables index_group;
run;

/*
 index_                             Cumulative    Cumulative
 group     Frequency     Percent     Frequency      Percent
 -----------------------------------------------------------
 HIGH            41       77.36            41        77.36
 LOW              9       16.98            50        94.34
 NONE             3        5.66            53       100.00
*/

/* Quick check of the "NONE" index_group:

Index SNPs without LD data (and candidate gene according to paper):
rs2611215	intergenic region (4q32.3)
rs35667974	IFIH1
rs689		INS
*/

/**************** CHECK LD BETWEEN T1D-ASSOC SNP and eQTL TESTED SNPs ****************/
/* Merge in the highest r2 value -- this will give us an idea about how well-linked 
 O-G credible SNPs are to the eQTL-tested SNPs */

proc sort data=eqtl.onengut_snps_ld_eqtl_snps;
   by onengut_index_snp descending r2;
run;

data onengut_index_snps_r2;
   set eqtl.onengut_snps_ld_eqtl_snps; 
   by onengut_index_snp;
   if first.onengut_index_snp then output;
run;

proc sort data=onengut_index_snp_eqtl_check;
   by onengut_index_snp;
proc sort data=onengut_index_snps_r2;
   by onengut_index_snp;
run;

data onengut_index_snp_highest_r2;
   merge onengut_index_snp_eqtl_check (in=in1) onengut_index_snps_r2 (in=in2);
   by onengut_index_snp;
   if in1 and in2 then output;
   else if in1 then do;
      r2=0; output; end;
   else output;
run;

/* How many index SNPs with r2<0.5, <0.7, <0.8? */

data onengut_index_snp_r2_flags;
   set onengut_index_snp_highest_r2;
   if r2 lt 0.5 then flag_r2_lt05=1; else flag_r2_lt05=0;
   if r2 lt 0.7 then flag_r2_lt07=1; else flag_r2_lt07=0;
   if r2 lt 0.8 then flag_r2_lt08=1; else flag_r2_lt08=0;
   if r2 lt 0.9 then flag_r2_lt09=1; else flag_r2_lt09=0;
run;

ods listing; ods html close;
proc freq data=onengut_index_snp_r2_flags;
   tables flag_r2_lt05 flag_r2_lt07 flag_r2_lt08 flag_r2_lt09;
run;

*41 index SNPs lt 0.9, 0.8, 0.7, 0.5;
*12 index SNPs lt 0.9, 0.8, 0.7, 0.5 (includes the 3 index SNPs without LD data);

data onengut_index_snp_low_ld;
   set onengut_index_snp_r2_flags;;
   if flag_r2_lt09=1;
   keep onengut_index_snp r2;
run;

*merge in likely candidate gene;

data index_snp2gene;
   set eqtl.onengut_supptable1;
   keep index_snp_rs genes;
   rename index_snp_rs=onengut_index_snp genes=gene_id;
run;

proc sort data=onengut_index_snp_low_ld;
   by onengut_index_snp;
proc sort data=index_snp2gene nodup;
   by onengut_index_snp gene_id;
run;

data index_snp_low_ld_w_gene;
   merge onengut_index_snp_low_ld (in=in1) index_snp2gene (in=in2);
   by onengut_index_snp;
   if in1 and in2;
run;

proc print data=index_snp_low_ld_w_gene;
run;

/* Index SNPs are:
         onengut_
  Obs    index_snp                r2    gene_id	candidate gene (from Onengut paper, Table 1)

    1    rs10277986        0.0494222
    2    rs1456988          0.242889
    3    rs2611215                 0
    4    rs34536443        0.0124088    TYK2	TYK2
    5    rs35667974                0		IFIH1
    6    rs35667974                0    IFIH1	IFIH1
    7    rs61839660        0.0615954    IL2RA	IL2RA
    8    rs6691977          0.209739
    9    rs689                     0    INSIGF2	INS
   10    rs72853903         0.183564    MIR4686	INS
   11    rs8056814          0.291015		BCAR1
   12    rs8056814          0.291015    CTRB1	BCAR1
   13    rs8056814          0.291015    CTRB2	BCAR1
   14    rs9585056          0.287178		GPR183


*/

/* Check MAFs for these index SNPs and the list of statistically-indistinguishable SNPs */

data index_snps_for_maf_check;
   set onengut_index_snp_low_ld;
   keep onengut_index_snp;
run;

data index_snp2credible_no_ld;
   set eqtl.onengut_snps_ld_no_eqtl_snps;
   *set observations were ld_test_snp_id is blank to their index snp;
   if ld_test_snp_id='' then ld_test_snp_id=onengut_index_snp;
   keep ld_test_snp_id onengut_index_snp;
run;

data index_snp2credible_ld;
   set eqtl.onengut_snps_ld_eqtl_snps;
   keep test_onengut_snp_id onengut_index_snp;
   rename test_onengut_snp_id=ld_test_snp_id;
run;

data index_snp2credible;
   set index_snp2credible_ld index_snp2credible_no_ld;
run;

proc sort data=index_snps_for_maf_check nodup;
   by onengut_index_snp;
proc sort data=index_snp2credible nodup;
   by onengut_index_snp;
run;

data credible_snps_for_maf_check no_credible_oops;
   merge index_snp2credible (in=in1) index_snps_for_maf_check (in=in2);
   by onengut_index_snp;
   if in1 and in2 then output credible_snps_for_maf_check;
   else if in2 then output no_credible_oops; *0 obs!;
run;

*407 SNPs to check minor allele freq;
* Get MAFs;

data snp_mafs;
   set eqtl.subset_maf_results;
   keep snp_id maf;
   rename snp_id=ld_test_snp_id;
run;

proc sort data=snp_mafs nodup;
   by ld_test_snp_id;
proc sort data=credible_snps_for_maf_check nodup;
   by ld_test_snp_id;
run;

data credible_snps_w_maf no_maf;
   merge snp_mafs (in=in1) credible_snps_for_maf_check (in=in2);
   by ld_test_snp_id;
   if in1 and in2 then output credible_snps_w_maf;
   else if in2 then output no_maf;
run;

*79 SNPs with MAFs;
/* 2 SNPs without MAFs: rs689 and rs34536443
Quick check with the ImmunoChip data: none appear in the genotyping data -- talk to Pat about this */

/* Put MAFs into bins to see how many are within certain MAF ranges:
<5%, 5-10%, 10-20%, 20-30%, 30-40%, 40-50% */

data credible_snps_maf_bins;
   set credible_snps_w_maf;
   length maf_bin $5.;
   if maf lt 0.05 then maf_bin='lt_05'; 
   else if maf lt 0.10 then maf_bin='05_10';
   else if maf lt 0.20 then maf_bin='10_20';
   else if maf lt 0.30 then maf_bin='20_30';
   else if maf lt 0.40 then maf_bin='30_40';
   else if maf le 0.50 then maf_bin='40_50';
   else maf_bin='gt_50';
run;

proc freq data=credible_snps_maf_bins;
    tables maf_bin;
run;

/* 
                                     Cumulative    Cumulative
 maf_bin    Frequency     Percent     Frequency      Percent
 ------------------------------------------------------------
 05_10             5        6.33             5         6.33
 10_20            11       13.92            16        20.25
 20_30            10       12.66            26        32.91
 30_40            27       34.18            53        67.09
 40_50             1        1.27            54        68.35
 lt_05            25       31.65            79       100.00
*/

data low_maf_indices high_maf_indices;
   set credible_snps_maf_bins;
   if maf_bin="lt_05" then output low_maf_indices;
   else output high_maf_indices;
   keep onengut_index_snp;
run;

proc sort data=low_maf_indices nodup;
   by onengut_index_snp;
proc sort data=high_maf_indices nodup;
   by onengut_index_snp;
run;

/* Check if these SNPs were tested for eQTLs */

data eqtl_snps;
   set eqtl.eqtl_results_summary_table;
   keep snp_id gene_id;
run;

data low_ld_snps_w_maf;
   set credible_snps_w_maf;
   rename ld_test_snp_id=snp_id;
run;

proc sort data=eqtl_snps;
   by snp_id;
proc sort data=low_ld_snps_w_maf;
   by snp_id;
run;

data low_ld_snps_w_eqtl;
   merge low_ld_snps_w_maf (in=in1) eqtl_snps (in=in2);
   by snp_id;
   if in1 and in2;
run;
* 0 SNPs were tested as eQTLs;
* Need to find why. Need to match SNP to gene;
* Can't do with SNPs <0.05 though.

/* Get list of SNPs to check */

data snps_for_gene;
   set low_ld_snps_w_maf;
   if maf ge 0.05;
run;

data snp_data_from_info;
  set eqtl.snp_data_w_info;
  keep snp_id gene_id;
run;

proc sort data=snp_data_from_info nodup;
   by snp_id gene_id;
proc sort data=snps_for_gene nodup;
   by snp_id;
run;


data snps_check_gene no_gene;
   merge snps_for_gene (in=in1) snp_data_from_info (in=in2);
   by snp_id;
   if in1 and in2 then output snps_check_gene;
   else if in1 then output no_gene;
run;

* Okay, these 54 SNPs are likely outside 5kb for a given gene;

/* Check what genes they're assigned to by the O-G paper */

data onengut_credible_snp2gene;
   set eqtl.onengut_supptable1;
   keep snp_id index_snp_rs cred_snp_rs genes;
run;

proc sort data=snps_for_gene;
   by snp_id;
proc sort data=onengut_credible_snp2gene;
   by snp_id;
run;

data snps_check_og_gene no_gene;
   merge snps_for_gene (in=in1) onengut_credible_snp2gene (in=in2);
   by snp_id;
   if in1 and in2 then output snps_check_og_gene;
   else if in1 then output no_gene;
run;

* Most SNPs are not assigned a gene according to Onengut paper;
* The few that have a gene are in either CTRB1, CTRB2, MIR4686;
* Check expression of these genes!;

data eqtl_genes;
   set eqtl.eqtl_results_summary_table;
   if index(gene_id, 'CTRB1') ge 1
   or index(gene_id, 'CTRB2') ge 1
   or index(gene_id, 'MIR4686') ge 1;
   keep gene_id;
run;

/* Okay, these genes not expressed! */

/* Make permenant sets of SNPs to keep for further analysis */


/* Summary:
3 O-G index SNPs have no LD data to the SNPs selected for eQTLs:
rs2611215 (4q32.3), rs35667974 (IFIH1) and rs689 (INS)
These index SNPs are therefore excluded from eQTL analysis

41 index SNPs in the "HIGH" LD group, 9 index SNPs in the "LOW" LD group



I then looked to see how many index SNP groups (this is the index SNP plus the statistically-indistinguishable SNPs) demonstrated evidence of high linkage with an eQTL-tested SNP:
41 index SNP groups had high (r2>0.9)
12 groups had low (r2<0.5): this includes the three index SNPs above.

Excluding the 3 with no LD, and removing duplicated index SNPs, this actually leaves us with 8 index SNP groups with low LD. These are:
rs10277986, rs1456988, rs34536443 (TYK2), rs61839660 (IL2RA), rs6691977, rs72853903 (MIR4686/INS), rs8056814 (BCAR1/CTRB1/CTRB2), and rs9585956 (GPR183)

Taking the SNPs from these low/no LD SNP groups gives 81 SNPs. I then checked the MAFs of these SNPs and found that 2 SNPs did not have MAF info from the available genotype data (rs689, rs34536443; both in INS). Of the remaining 79 SNPs with MAF info, 25 had a MAF <5%. These were from the index groups: rs34536443, rs10277986, rs61839660, rs35667974. The 54 SNPs with MAFs >5% were from the index groups: rs1456988, rs2611215, rs6691977, rs72853903, rs8056814, rs9585056.

I then checked to see if the low/no LD SNPs were tested as eQTLs. None were.
I then checked the low/no LD SNPs with MAFs >5% to see if any of these were within the 5kb region for tested genes. None were extracted, suggesting that these SNPs fall outside the 5kb boundary set for each gene.

I then took the gene assignments for these 54 SNPs from the O-G paper and checked if these genes were expressed in our data. Most SNPs were not assigned a gene, but the few that were had the gene assignments: CTRB1, CTRB2, MIR4686. Checking the eQTL summary table suggested that none of these genes were expressed.

*/

