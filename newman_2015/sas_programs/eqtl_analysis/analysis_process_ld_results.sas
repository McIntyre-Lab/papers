/* Import LD data and process */

/* Set libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Import LD results - changing column headers to be more descriptive */

   data WORK.SUBSET_LD_RESULTS    ;
   %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
   infile '/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eurfam_subset_ld.csv'
   delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
      informat chr_ic_snp $4. ;
      informat pos_ic_snp best32. ;
      informat snp_id $15. ;
      informat chr_onengut_snp $4. ;
      informat pos_onengut_snp best32. ;
      informat test_onengut_snp_id $15. ;
      informat r2 best32. ;
      format chr_ic_snp $4. ;
      format pos_ic_snp best12. ;
      format snp_id $15. ;
      format chr_onengut_snp $4. ;
      format pos_onengut_snp  best12. ;
      format test_onengut_snp_id $15. ;
      format r2 best12. ;
   input
               chr_ic_snp $
               pos_ic_snp
               snp_id $
               chr_onengut_snp $
               pos_onengut_snp 
               test_onengut_snp_id $
               r2
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

/* Import minor allele frequency data - need this for later */

    data WORK.SUBSET_MAF_RESULTS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '/home/jrbnewman/concannon/eqtl_analysis/pipeline_output/eurfam_subset_maf.csv'
delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat chr $4. ;
       informat snp_id $15. ;
       informat minor_allele $1. ;
       informat major_allele $1. ;
       informat maf best32. ;
       informat num_alleles_obs best32. ;
       format chr $4. ;
       format snp_id $15. ;
       format minor_allele $1. ;
       format major_allele $1. ;
       format maf best12. ;
       format num_alleles_obs best12. ;
    input
                chr $
                snp_id $
                minor_allele $
                major_allele $
                maf
                num_alleles_obs
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;


/* Make permenant */

data eqtl.subset_ld_results;
   set subset_ld_results;
run;

data eqtl.subset_maf_results;
   set subset_maf_results;
run;

/* Check that the test SNPs (snp_id) are only the eQTL-tested SNPs */

data ld_snp_check;
   set eqtl.subset_ld_results;
   keep snp_id;
run;

data eqtl_tested_snps;
   set eqtl.eqtl_results_summary_table;
   keep snp_id;
run;

proc sort data=ld_snp_check nodup;
   by snp_id;
proc sort data=eqtl_tested_snps nodup;
   by snp_id;
run;

data ld_eqtl_test_snp_check no_ld not_eqtl_snp;
   merge eqtl_tested_snps (in=in1) ld_snp_check (in=in2);
   by snp_id;
   if in1 and in2 then output ld_eqtl_test_snp_check;
   else if in1 then output no_ld;
   else output not_eqtl_snp;
run;

*11368 SNPs in eQTL-tested list, 11368 SNPs in LD list;
*11368 SNPs with LD data, 0 without LD data, 0 are not eQTL-tested SNPs!;


/* Need to merge in the index SNP with the LD SNPs */

/* How many SNPs have multiple index SNPs?? */

data index_snp_count_check;
   set eqtl.supptable1_snp_list;
   keep snp_id index_snp_rs;
run;

proc sort data=index_snp_count_check nodup;
   by snp_id index_snp_rs;
run;

proc freq data=index_snp_count_check noprint;
    tables snp_id / out=index_snp_count;
run;

data index_snp_count_multi;
   set index_snp_count;
   if count gt 1;
run;

*4 SNPs with 2 index SNPs, but at same genes. Will cat index SNPs together ;
*max=2 SNPs per extracted SNP;

data cat_index_snps; 
  array index_snp[2] $ 15.;
  retain index_snp1-index_snp2;
  set index_snp_count_check;
  by snp_id;
  if first.snp_id then do;
     call missing(of index_snp1-index_snp2);
     records = 0;
  end;
  records + 1;
  index_snp[records]=index_snp_rs;
  if last.snp_id then output;
run;

  *clean up the output file;

data cat_index_snps_2;
  set cat_index_snps; 
  length onengut_index_snp $ 31.;
         onengut_index_snp= catx("|", OF index_snp1-index_snp2);
  drop index_snp1-index_snp2 index_snp_rs records;
  run;

/* Convert O-G SNP IDs to ImmunoChip SNP IDs */

data onengut2immunochip_id;
   set eqtl.ic_to_onengut_snp_index;
   keep snp_id onengut_snp_id;
run;

data onengut_assoc_snps;
   set cat_index_snps_2;
   rename snp_id=onengut_snp_id;
run;

proc sort data=onengut2immunochip_id;
   by onengut_snp_id;
proc sort data=onengut_assoc_snps;
   by onengut_snp_id;
run;

data og_assoc_snps_ic_ids no_assoc no_ic_id_oops;
   merge onengut2immunochip_id (in=in1) onengut_assoc_snps (in=in2);
   by onengut_snp_id;
   if in1 and in2 then output og_assoc_snps_ic_ids;
   else if in1 then output no_assoc;
   else output no_ic_id_oops;
run;

/* Merge Onengut SNPs to LD SNPs */

data subset_ld_results;
   set eqtl.subset_ld_results;
run;

data onengut_w_ic_ids;
   set og_assoc_snps_ic_ids;
   rename snp_id=test_onengut_snp_id;
run;

proc sort data=subset_ld_results;
    by test_onengut_snp_id snp_id;
proc sort data=onengut_w_ic_ids nodup;
    by test_onengut_snp_id;
run;

data onengut_ic_snp_w_eqtl_snp no_eqtl no_onengut;
   merge onengut_w_ic_ids (in=in1) subset_ld_results (in=in2);
   by test_onengut_snp_id;
   if in1 and in2 then output onengut_ic_snp_w_eqtl_snp;
   else if in1 then output no_eqtl;
   else output no_onengut;
run;

*2027 Onengut_2015 SNPs, 207897 LD results;
* 7958 LD results with an Onengut SNP;
* 929 Onengut SNPs do not have an LD SNP;
* 199939 LD results without an Onengut SNP;

/* Check that remaining LD results are between eQTL-tested SNPs */

data eqtl_tested_snps;
   set eqtl.eqtl_results_summary_table;
   keep snp_id;
   rename snp_id=test_onengut_snp_id;
run;

data remaining_ld_results;
  set no_onengut;
  drop onengut_snp_id onengut_index_snp;
run;

proc sort data=remaining_ld_results;
    by test_onengut_snp_id snp_id;
proc sort data=eqtl_tested_snps nodup;
    by test_onengut_snp_id;
run;

data ld_Results_eqtl ld_results_non_eqtl_oops no_ld;
   merge eqtl_tested_snps (in=in1) remaining_ld_results (in=in2);
   by test_onengut_snp_id;
   if in1 and in2 then output ld_Results_eqtl;
   else if in1 then  output no_ld;
   else output ld_results_non_eqtl_oops; *0 obs, yay!;
run;

*11368 eQTL-tested SNPs, 199939 SNPs remaining from LD results;
* 199939 LD SNPs were tested as eQTLs, 0 remaining LD results from non-eQTL-tested SNPs, 145 eQTL-SNPs with LD results;
* all remaining LD results are between eQTL SNPs;

/* Make permenant */

data eqtl.onengut_snps_ld_eqtl_snps;
   set onengut_ic_snp_w_eqtl_snp;
run;

data eqtl.onengut_snps_ld_no_eqtl_snps;
   set no_eqtl;
   keep test_onengut_snp_id onengut_snp_id onengut_index_snp;
   rename test_onengut_snp_id=ld_test_snp_id;
run;


/* TO DO:

1. check if any Onengut index SNP does not have a match to an eQTL SNP
	- This should give the same result as the previous matchin
2. make the list of credible snps: high (r2>80%), moderate (r2>70%), modest (r2>70%)
3. Merge in the eQTL results
4. Count number of Index SNPs (ie, regions) with an eQTL tested
5. Count number of Index SNPs (ie, regions) with an eQTL significant
6. Make a table for exporting */


