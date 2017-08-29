/*find:
eqtl_summary_by_index_snp.csv (for Table 5)
make_summary_table_for_index_snp.sas

Supplemental_Table_1.csv (for table S3 -- eqtL results all)
analysis_eqtl_update_diabetes_flag.sas
*/

/* update main table, source from this for index SNP summaries 

   For each, I will have two magnitudes: HTZ vs HMZ1, HMZ2 vs HMZ1
   Then for each, take the HMZ2/1 where possible, otherwith HTZ/HMZ1

*/

libname eqtl '!PATCON/eqtl_analysis/sas_data';
libname con '!PATCON/sas_data';
libname eqtltmp '/mnt/store/eqtl_sandbox';

/* For each SNP, get the counts per genotype, as well as the major and minor alleles */

data gt_cnts;
  set eqtl.genotype_data;
  if subject_id="M080" then delete;
  if subject_id="M068" then delete;
  keep snp_id genotype subject_id;
run;

proc sort data=gt_cnts nodup;
   by snp_id genotype subject_id;
proc freq data=gt_cnts noprint;
   by snp_id;
   tables genotype / out=gt_cnts_by_snp;
run;

data gt_cnts_by_snp2;
  length genotype2 $20.;
  set gt_cnts_by_snp;
  if genotype=0 then genotype2="num_hmz1";
  if genotype=1 then genotype2="num_htz";
  if genotype=2 then genotype2="num_hmz2";
  if genotype=. then delete;
run;

proc transpose data=gt_cnts_by_snp2 out=gt_cnts_sbys;
   by snp_id;
   var count;
   id genotype2;
run;

data gt_cnts_sbys2;
   set gt_cnts_sbys;
   array change _numeric_;
      do over change;
      if change=. then change=0;
      end;
   drop _NAME_ _LABEL_;
run;

* Alleles;

data alleles;
  set eqtl.snp_data_w_info;
  keep snp_id ref_allele alt_allele;
  rename ref_allele=allele1 alt_allele=allele2;
run;

proc sort data=alleles nodup;
   by snp_id;
proc sort data=gt_cnts_sbys2;
   by snp_id;
run;

data gt_count_w_info;
  merge alleles (in=in1) gt_cnts_sbys2 (in=in2);
  by snp_id;
  if in1 and in2;
run;

  
/* First, transpose data so that I can make estimates +/i SD, and calculate diffs (for magnitude) */

data lsmeans;
   set eqtltmp.eqtl_lsmeans_results;
   keep feature_id snp_id cell_type genotype estimate stderr;
run;

proc sort data=lsmeans nodup;
   by feature_id snp_id cell_type genotype estimate stderr;
proc transpose data=lsmeans out=lsmeans_sbys_est;
   by feature_id snp_id cell_type;
   var estimate;
   id genotype;
run;
proc transpose data=lsmeans out=lsmeans_sbys_stderr;
   by feature_id snp_id cell_type;
   var stderr;
   id genotype;
run;

data lsmeans_sbys_est2;
   set lsmeans_sbys_est;
   drop _NAME_;
   rename _0=hmz1_estimate _1=htz_estimate _2=hmz2_estimate;
run;
 

data lsmeans_sbys_stderr2;
   set lsmeans_sbys_stderr;
   drop _NAME_ _LABEL_;
   rename _0=hmz1_se _1=htz_se _2=hmz2_se;
run;


proc sort data=lsmeans_sbys_est2;
   by feature_id snp_id cell_type;
proc sort data=lsmeans_sbys_stderr2;
   by feature_id snp_id cell_type;
run;


data lsmeans_sbys;
  merge lsmeans_sbys_est2 (in=in1) lsmeans_sbys_stderr2 (in=in2);
  by feature_id snp_id cell_type;
  if in1 and in2;
run;

data calc_magnitude;
  set lsmeans_sbys;
  if hmz1_estimate > htz_estimate then magnitude_htz_hmz1=-hmz1_estimate/htz_estimate; else magnitude_htz_hmz1=htz_estimate/hmz1_estimate;
  if hmz1_estimate > hmz2_estimate then magnitude_hmz2_hmz1=-hmz1_estimate/hmz2_estimate; else magnitude_hmz2_hmz1=hmz2_estimate/hmz1_estimate;
run;

data estimate_range;
  format hmz1_estimate f8.3;
  format htz_estimate f8.3;
  format hmz1_estimate f8.3;
  format hmz1_se f8.3;
  format htz_se f8.3;
  format hmz1_se f8.3;
  set calc_magnitude;
  if hmz1_estimate ne . then do;
      hmz1_estimate_range=catt(hmz1_estimate, ' ± ', hmz1_se);
      htz_estimate_range=catt(htz_estimate, ' ± ', htz_se);
      hmz2_estimate_range=catt(hmz2_estimate, ' ± ', hmz2_se);
      end;
  else do;
      hmz1_estimate_range="n.d.";
      htz_estimate_range="n.d.";
      hmz2_estimate_range="n.d.";
  end;
  keep feature_id snp_id cell_type hmz1_estimate_range htz_estimate_range hmz2_estimate_range
       magnitude_htz_hmz1 magnitude_hmz2_hmz1;
  rename hmz1_estimate_range=hmz1_estimate
         htz_estimate_range=htz_estimate
         hmz2_estimate_range=hmz2_estimate;
run;

data cd4 cd8 cd19;
   set estimate_range;
   if cell_type="CD4" then output cd4;
   if cell_type="CD8" then output cd8;
   if cell_type="CD19" then output cd19;
   drop cell_type;
run;

data cd4_2;
  set cd4;
  rename hmz1_estimate=hmz1_estimate_CD4
         htz_estimate=htz_estimate_CD4
         hmz2_estimate=hmz2_estimate_CD4
         magnitude_htz_hmz1=magnitude_htz_hmz1_CD4
         magnitude_hmz2_hmz1=magnitude_hmz2_hmz1_CD4;
run;

data cd8_2;
  set cd8;
  rename hmz1_estimate=hmz1_estimate_CD8
         htz_estimate=htz_estimate_CD8
         hmz2_estimate=hmz2_estimate_CD8
         magnitude_htz_hmz1=magnitude_htz_hmz1_CD8
         magnitude_hmz2_hmz1=magnitude_hmz2_hmz1_CD8;
run;

data cd19_2;
  set cd19;
  rename hmz1_estimate=hmz1_estimate_CD19
         htz_estimate=htz_estimate_CD19
         hmz2_estimate=hmz2_estimate_CD19
         magnitude_htz_hmz1=magnitude_htz_hmz1_CD19
         magnitude_hmz2_hmz1=magnitude_hmz2_hmz1_CD19;
run;

proc sort data=cd4_2;
  by feature_id snp_id;
proc sort data=cd8_2;
  by feature_id snp_id;
proc sort data=cd19_2;
  by feature_id snp_id;
run;

data results_summary;
  set eqtl.results_summary_table_w_means_v2;
run;


data results_summary_cred;
  set eqtl.results_w_onengut_cred_means_v2;
run;

proc sort data=results_summary;
   by feature_id snp_id;
proc sort data=results_summary_cred;
  by feature_id snp_id;
run;

data results_summary_w_est;
  merge results_summary (in=in1) cd4_2 cd8_2 cd19_2;
  by feature_id snp_id;
  if in1;
run;

data results_summary_cred_w_est;
  merge results_summary_cred (in=in1) cd4_2 cd8_2 cd19_2;
  by feature_id snp_id;
  if in1;
run;

proc sort data=results_summary_w_est;
   by snp_id;
proc sort data=results_summary_cred_w_est;
   by snp_id;
proc sort data=gt_count_w_info;
   by snp_id;
run;


data results_summary_w_est2;
  merge gt_count_w_info (in=in1) results_summary_w_est (in=in2);
  by snp_id;
  if in1 and in2;
run;

data results_summary_cred_w_est2;
  merge gt_count_w_info (in=in1) results_summary_cred_w_est (in=in2);
  by snp_id;
  if in1 and in2;
run;

data eqtl.results_summary_table_w_means_v4;
   retain gene_id flag_multigene flag_autoimmune_gene flag_diabetes_gene feature_type feature_id snp_id allele1 allele2 num_hmz1 num_htz num_hmz2;
   set results_summary_w_est2;
run;

data eqtl.results_w_onengut_cred_means_v4;
   retain gene_id snp_id allele1 allele2 num_hmz1 num_htz num_hmz2;
   set results_summary_cred_w_est2;
run;

proc export data=eqtl.results_summary_table_w_means_v3
     outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/all_eqtl_results1.csv' dbms=csv replace;
run;

proc export data=eqtl.results_w_onengut_cred_means_v3
     outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/all_eqtl_results_w_index_snps1.csv' dbms=csv replace;
run;




proc export data=eqtl.results_summary_table_w_means_v4
     outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/all_eqtl_results2.csv' dbms=csv replace;
run;

proc export data=eqtl.results_w_onengut_cred_means_v4
     outfile='/home/jrbnewman/McLab/jrbnewman/manuscripts/Newman_T1D_splicing/reviewer_responses/all_eqtl_results_w_index_snps2.csv' dbms=csv replace;
run;






