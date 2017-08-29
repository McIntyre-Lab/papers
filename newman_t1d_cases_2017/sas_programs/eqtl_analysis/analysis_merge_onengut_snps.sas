libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Now I need to merge the list of tested SNPs with the actual results! */

/* First check that each snp entry is unique. Otherwise, merges will be bad!! */

data onengut_snps_w_tests;
   set eqtl.onengut_snps_w_tests;
   keep ucsc_snp_id snp_id tested_snp_id;
run;

proc sort data=onengut_snps_w_tests;
   by snp_id;
run;
* check of data looks like some SNPs could be assigned to multiple tested SNPs;
* will concatenate the Onengut SNP IDs per tested SNP ID for ease;


data onengut_snps_to_cat;
   set onengut_snps_w_tests;
   keep snp_id tested_snp_id;
run;

proc sort data=onengut_snps_to_cat;
   by tested_snp_id snp_id;
run;


/*** cat snp_ids ***/

/* get counts first */
proc freq noprint data=onengut_snps_to_cat;
   tables tested_snp_id / out=snp_count;
run;

proc sort data=snp_count;
  by descending count;
run;
*max=52 SNPs per extracted SNP;

data cat_snps; 
  array snp[52] $ 15.;
  retain snp1-snp52;
  set onengut_snps_to_cat;
  by tested_snp_id;
  if first.tested_snp_id then do;
     call missing(of snp1-snp52);
     records = 0;
  end;
  records + 1;
  snp[records]=snp_id;
  if last.tested_snp_id then output;
run;

  *clean up the output file;

data cat_snps_2;
  set cat_snps;
  length onengut_snp_ids $ 832.;
  rename records= num_cred_snps;
         onengut_snp_ids= catx("|", OF snp1-snp52);
  drop snp1-snp52 snp_id;
  rename tested_snp_id=snp_id;
  run;


/* Merge in with the signficant results */

data eqtl_results;
   set eqtl.eqtl_results_summary_table;
run;

proc sort data=eqtl_results;
   by snp_id;
proc sort data=cat_snps_2;
   by snp_id;
run;

data eqtl_results_w_assoc_snps no_assoc no_results;
   merge eqtl_results (in=in1) cat_snps_2 (in=in2);
   by snp_id;
   if in1 and in2 then output eqtl_results_w_assoc_snps;
   else if in1 then output no_assoc;
   else output no_results;
run;



/* Check of the ones with assoc, how many have eqtls */

data eqtl_results_assoc_snps;
   set eqtl_results_w_assoc_snps;
   if flag_eqtl_sig=1;
run;




rs1893592

   else if in1 then do;
        num_cred_snps=0;
        onengut_snp_ids='';
        output eqtl_results_w_assoc_snps; end;












* Need a Onengut SNP to Region table too, for checking genes etc.;
data eqtl.onengut_snp2gene;
   set eqtl.onengut_supptable1;
   keep snp_id index_snp_rs cred_snp_rs genes;
run;
