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

data counts;
  set event.t1d_flag_mei_cell_sbys;
run;

proc sort data=genes_for_plots;
  by gene_id;
proc sort data=counts;
  by gene_id;
run;

data counts2;
  merge genes_for_plots (in=in1) counts (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=counts2;
  by gene_id transcript_id;
proc means data=counts2 noprint;
  by gene_id transcript_id;
  var TPM_CD19 TPM_CD4 TPM_CD8;
  output out=mean_counts mean=;
run;

data export_counts;
  set mean_counts;
  drop _TYPE_ _FREQ_;
run;

proc export data=export_counts outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_mean_tpm_diff_mei.csv"
   dbms=csv replace;
run;

*output gene list;

data genes;
   set export_counts;
   keep gene_id;
run;

proc sort data=genes nodup;
  by gene_id;
run;

proc export data=genes outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_gene_list_diff_mei.csv"
   dbms=csv replace; putnames=no;
run;


