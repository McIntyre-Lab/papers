ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* For T1D genes with splicing differences as per the GLM model test, export data and make some plots
   I want both the "heatmap" plot (with MEE counts) and the fragment plot */

data t1d_genes;
  set con.immunobase_gene_flags;
  where flag_immunobase_diabetes_gene=1;
  keep gene_id;
run;

data sig_genes;
  set event.t1d_rank_glm_results_by_gene;
  where flag_mee_p05=1;
  keep gene_id;
run;

data fus2gene;
   set hg19.hg19_av_fusion_si_info_unique;
   keep fusion_id gene_id fusion_start;
run;

proc sort data=fus2gene nodup;
  by gene_id;
proc sort data=t1d_genes nodup;
  by gene_id;
proc sort data=sig_genes nodup;
  by gene_id;
run;

data fus2keep;
  merge fus2gene (in=in1) t1d_genes (in=in2) sig_genes (in=in3);
  by gene_id;
  if in1 and in2 and in3;
run;


data fus_counts;
  set con.fusion_q3_norm_data_all;
run;

proc sort data=fus2keep nodup;
  by fusion_id gene_id;
proc sort data=fus_counts;
  by fusion_id;
run;

data fus_counts2keep;
  merge fus2keep (in=in1) fus_counts (in=in2);
  by fusion_id;
  if in1 and in2;
run;

data design;
  set con.design_by_subject_new;
  if subject_id in ("M075","M048","M080") then delete;
  keep name subject_id cell_Type;
run;

proc sort data=fus_counts2keep;
  by name;
proc sort data=design;
  by name;
run;

data fus2gene_counts_w_key;
  merge design (in=in1) fus_counts2keep (in=in2);
  by name;
  if in1 and in2;
run;


proc sort data=fus2gene_counts_w_key;
   by gene_id subject_id fusion_start fusion_id cell_type;
proc transpose data=fus2gene_counts_w_key out=counts_sbys;
   by gene_id subject_id fusion_start  fusion_id ;
   var log_q3_q3_apn;
   id cell_type;
run;

data counts;
  set counts_sbys;
  drop _NAME_;
  rename cd4=logapn_cd4 cd8=logapn_cd8 cd19=logapn_cd19;
run;

proc sort data=counts;
  by gene_id fusion_start fusion_id;
proc means data=counts noprint;
  by gene_id fusion_start fusion_id;
  var logapn_CD19 logapn_CD4 logapn_CD8;
  output out=mean_counts mean=;
run;

data mean_counts2;
  set mean_counts;
  drop _TYPE_ _FREQ_;
run;

/* Count times where a transcript is the MEI */

proc sort data=fus2gene_counts_w_key;
   by gene_id cell_type subject_id descending log_q3_q3_apn;
run;

data flag_mee;
  set fus2gene_counts_w_key;
  by gene_id cell_type subject_id;
  if first.subject_id then flag_mee=1;
  else flag_mee=0;
run;

proc sort data=flag_mee;
  by gene_id fusion_id cell_type ;
proc means data=flag_mee noprint;
  by gene_id fusion_id cell_type ;
   var flag_mee;
  output out=sum_flag_mee sum=;
run;

proc transpose data=sum_flag_mee out=mee_sbys;
  by gene_id fusion_id;
  id cell_type;
  var flag_mee;
run;

data mee_sbys2;
  set mee_sbys;
  drop _NAME_;
  rename CD19=flag_mee_cd19 CD4=flag_mee_cd4 CD8=flag_mee_cd8;
run;

proc sort data=mee_sbys2;
   by fusion_id;
proc sort data=mean_counts2;
   by fusion_id;
run;

data export_counts_mee;
  merge mean_counts2 (in=in1) mee_sbys2 (in=in2);
  by fusion_id;
  if in1 and in2 then output;
  else if in1 then do;
     flag_mee_cd19=0;
     flag_mee_cd4=0;
     flag_mee_cd8=0;
     output; end;
run;

proc sort data=export_counts_mee;
  by gene_id fusion_start;
run;

data export_counts_mee2;
  set export_counts_mee;
  drop fusion_start;
run;

data gene_list;
  set export_counts_mee2;
  keep gene_id;
run;

proc sort data=gene_list nodup;
  by gene_id;
run;


proc export data=export_counts_mee2 outfile="!MCLAB/event_analysis/analysis_output/t1d_genes_MEE_glm_test_for_heatmaps.csv"
   dbms=csv replace;
run;

proc export data=gene_list outfile="!MCLAB/event_analysis/analysis_output/t1d_genes_MEE_glm_test_gene_list.csv"
   dbms=csv replace;
run;


