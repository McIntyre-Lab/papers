/* Set libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

/* Now we want to make a summary table that looks like this:
gene_id feature_id snp_id [eqtl_results] onengut_index_snp credible_snp_id r2

Since this can and will cause some issues with a many-to-many merge, I will have to concatentate the credible snp IDs together first then merge in the LD info
*/


/* Get the list of O-G SNPs to eQTL SNPs */
data onengut2eqtl_snps;
  set eqtl.onengut_snps_ld_eqtl_snps;
  keep test_onengut_snp_id snp_id;
run;

* check for any duplicates;
proc sort data=onengut2eqtl_snps nodup;
  by snp_id test_onengut_snp_id;
run;
* no dups!*;

/* Need to cat together test_onengut_snp_id for each snp_id. How many O-G SNPs per eQTL SNP? */

proc freq data=onengut2eqtl_snps noprint;
   tables snp_id /out=og_snp_count;
run;

proc sort data=og_snp_count;
   by descending count;
run;

*19 SNPs max;

* Cat together;

data cat_og_snps; 
  array og_snp[19] $ 15.;
  retain og_snp1-og_snp19;
  set onengut2eqtl_snps;
  by snp_id;
  if first.snp_id then do;
     call missing(of og_snp1-og_snp19);
     records = 0;
  end;
  records + 1;
  og_snp[records]=test_onengut_snp_id;
  if last.snp_id then output;
run;

  *clean up the output file;

data cat_og_snps_2;
  set cat_og_snps; 
  length test_og_snp_id $ 500.;
         test_og_snp_id= catx("|", OF og_snp1-og_snp19);
  keep snp_id test_og_snp_id;
  run;

/* Merge in eQTL summary */

data eqtl_summaries;
   set eqtl.eqtl_results_summary_table;
run;

proc sort data=eqtl.eqtl_results_summary_table out=eqtl_summaries;
   by snp_id;
proc sort data=cat_og_snps_2;
   by snp_id;
run;

data eqtl_results_t1d_snps not_t1d no_eqtl; 
   merge cat_og_snps_2 (in=in1) eqtl_summaries (in=in2);
   by snp_id;
   if in1 and in2 then output eqtl_results_t1d_snps;
   else if in1 then output no_eqtl; *0 obs, yay!;
   else output not_t1d;
run;

/* Make this permenant for now. May want to look at it again later */

data eqtl.eqtl_results_t1d_snps;
   set eqtl_results_t1d_snps;
run;


/* Un-concatenate O-G SNPs */

data eqtl_results_t1d_snps_stack;
   length test_onengut_snp_id $15.; 
   set eqtl_results_t1d_snps;
   do i=1 by 1 while(scan(test_og_snp_id,i,'|') ^=' ');
      test_onengut_snp_id=scan(test_og_snp_id,i,'|');
      output;
   end;
   drop i test_og_snp_id;
run;


/* Merge in LD results */
data onengut2eqtl_snps_ld;
  set eqtl.onengut_snps_ld_eqtl_snps;
  keep test_onengut_snp_id onengut_snp_id onengut_index_snp snp_id r2;
run;

proc sort data=eqtl_results_t1d_snps_stack;
   by snp_id test_onengut_snp_id;
proc sort data=onengut2eqtl_snps_ld;
   by snp_id test_onengut_snp_id;
run;


data eqtl_results_to_onengut_index no_onengut no_eqtl;
   merge eqtl_results_t1d_snps_stack (in=in1) onengut2eqtl_snps_ld (in=in2);
   by snp_id test_onengut_snp_id;
   if in1 and in2 then output eqtl_results_to_onengut_index;
   else if in1 then output no_onengut;
   else output no_eqtl;
run;

* rearrange variables and make permenant;
data eqtl_results_to_og_index;
   retain gene_id flag_diabetes_gene feature_id feature_type flag_multigene
   snp_id onengut_snp_id onengut_index_snp r2
   flag_eqtl_sig CD4_FDR_P CD8_FDR_P CD19_FDR_P;
   set eqtl_results_to_onengut_index;
   rename r2=onengut_snp2eqtl_snp_r2;
   drop flag_autoimmune_gene test_onengut_snp_id;
run;


data eqtl.eqtl_results_w_onengut_index;
   set eqtl_results_to_og_index;
run;


