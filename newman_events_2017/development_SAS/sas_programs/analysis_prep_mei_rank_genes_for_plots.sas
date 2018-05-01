/* Get genes that have very different MEIs based on the rank test
   for making cell heatmaps */
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';


data genes_for_plots;
   set event.t1d_case_only_kw_test;
   if flag_kw_cd4cd8_p05=1 and flag_kw_cd4cd19_p05=0 and flag_kw_cd8cd19_p05=0 then group_compare=1;
   if flag_kw_cd4cd8_p05=0 and flag_kw_cd4cd19_p05=1 and flag_kw_cd8cd19_p05=0 then group_compare=1;
   if flag_kw_cd4cd8_p05=0 and flag_kw_cd4cd19_p05=0 and flag_kw_cd8cd19_p05=1 then group_compare=1;
   if flag_kw_cd4cd8_p05=1 and flag_kw_cd4cd19_p05=1 and flag_kw_cd8cd19_p05=0 then group_compare=2;
   if flag_kw_cd4cd8_p05=1 and flag_kw_cd4cd19_p05=0 and flag_kw_cd8cd19_p05=1 then group_compare=2;
   if flag_kw_cd4cd8_p05=0 and flag_kw_cd4cd19_p05=1 and flag_kw_cd8cd19_p05=1 then group_compare=2;
   if flag_kw_cd4cd8_p05=1 and flag_kw_cd4cd19_p05=1 and flag_kw_cd8cd19_p05=1 then group_compare=3;
   keep gene_id group_compare;
run;

/*
65 genes different in only one pairwise comparison
206 genes different in two pairwise comparisons
54 genes different in all three pairwise comparisons

*/

data xs_counts;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
run;

data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  keep transcript_id gene_id;
run;

proc sort data=xs_counts;
  by transcript_id;
proc sort data=xs2gene;
  by transcript_id;
run;

data xs2gene_counts;
  merge xs2gene (in=in1) xs_counts (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data design;
  set con.design_by_subject_new;
  if subject_id in ("M075","M048","M080") then delete;
  keep library subject_id cell_Type;
run;

proc sort data=xs2gene_counts;
  by library;
proc sort data=design;
  by library;
run;

data xs2gene_counts_w_key;
  merge design (in=in1) xs2gene_counts (in=in2);
  by library;
  if in1 and in2;
run;

proc sort data=xs2gene_counts_w_key;
   by gene_id subject_id transcript_id cell_type;
proc transpose data=xs2gene_counts_w_key out=counts_sbys;
   by gene_id subject_id transcript_id ;
   var tpm;
   id cell_type;
run;

data counts;
  set counts_sbys;
  drop _NAME_;
  rename cd4=tpm_cd4 cd8=tpm_cd8 cd19=tpm_cd19;
run;

proc sort data=genes_for_plots;
  by gene_id;
proc sort data=counts;
  by gene_id;
run;

data counts2 nocount;
  merge genes_for_plots (in=in1) counts (in=in2);
  by gene_id;
  if in1 and in2 then output counts2;
  else if in1 then output nocount;
run;

proc sort data=counts2;
  by gene_id transcript_id;
proc means data=counts2 noprint;
  by gene_id transcript_id;
  var TPM_CD19 TPM_CD4 TPM_CD8;
  output out=mean_counts mean=;
run;

data mean_counts2;
  set mean_counts;
  drop _TYPE_ _FREQ_;
run;

/* Count times where a transcript is the MEI */

data mei;
   set event.t1d_flag_mei_cell_sbys;
   keep transcript_id subject_id flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
run;

proc sort data=mei;
   by transcript_id;
proc means data=mei noprint;
   by transcript_id;
   var flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
   output out=mei_count_by_xs sum=;
run;

data mei_count_by_xs2;
   set mei_count_by_xs;
   drop _TYPE_ _FREQ_;
run;

proc sort data=mei_count_by_xs2;
   by transcript_id;
proc sort data=mean_counts2;
   by transcript_id;
run;

data export_counts;
  merge mean_counts2 (in=in1) mei_count_by_xs2 (in=in2);
  by transcript_id;
  if in1 and in2 then output;
  else if in1 then do;
     flag_mei_cd19=0;
     flag_mei_cd4=0;
     flag_mei_cd8=0;
     output; end;
run;

proc export data=export_counts outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_mean_tpm_diff_mei_by_ranktest.csv"
   dbms=csv replace;
run;

*output gene list;

data genes1 genes2 genes3;
   set genes_for_plots;
   if group_compare=1 then output genes1;
   else if group_compare=2 then output genes2;
   else if group_compare=3 then output genes3;
   keep gene_id;
run;

proc sort data=genes1 nodup;
  by gene_id;
proc sort data=genes2 nodup;
  by gene_id;
proc sort data=genes3 nodup;
  by gene_id;
run;

proc export data=genes1 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_gene_list_ranktest_1comparison.csv"
   dbms=csv replace; putnames=no;
run;

proc export data=genes2 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_gene_list_ranktest_2comparisons.csv"
   dbms=csv replace; putnames=no;
run;

proc export data=genes3 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_gene_list_ranktest_3comparisons.csv"
   dbms=csv replace; putnames=no;
run;



