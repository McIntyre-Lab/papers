ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

/* Rank transcripts by TPM and perform K-W test pairwise by gene (CD4 vs CD8, CD4 vs CD19, CD8 vs CD19) */

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
   by subject_id cell_type gene_id descending tpm;
run;

proc rank data=xs2gene_counts_w_key out=xs_rank_by_gene_subj descending ties=low;
   by subject_id cell_type gene_id;
   var tpm ;
   ranks tpm_rank;
run;

proc sort data=xs_rank_by_gene_subj;
  by gene_id cell_type;
run;

proc npar1way data=xs_rank_by_gene_subj noprint ;
  by gene_id;
  class cell_type;
  var tpm_rank;
  output out=cd4cd8cd19_kw_xs_rank WILCOXON;
run;

proc npar1way data=xs_rank_by_gene_subj noprint ;
  where cell_Type ne "CD19";
  by gene_id;
  class cell_type;
  var tpm_rank;
  output out=cd4cd8_kw_xs_rank WILCOXON;
run;


proc npar1way data=xs_rank_by_gene_subj noprint ;
  where cell_Type ne "CD8";
  by gene_id;
  class cell_type;
  var tpm_rank;
  output out=cd4cd19_kw_xs_rank WILCOXON;
run;

proc npar1way data=xs_rank_by_gene_subj noprint ;
  where cell_Type ne "CD4";
  by gene_id;
  class cell_type;
  var tpm_rank;
  output out=cd8cd19_kw_xs_rank WILCOXON;
run;

data cd4cd8cd19_kw_2;
  set cd4cd8cd19_kw_xs_rank;
  KW_stat_all=_KW_;
  KW_stat_all_P=P_KW;
  if P_KW = . then flag_kw_all_p05=.;
  else if P_KW < 0.05 then flag_kw_all_p05=1;
  else flag_kw_all_p05=0;
  keep gene_id KW_stat_all KW_stat_all_P flag_kw_all_p05;
run;

data cd4cd8_kw_2;
  set cd4cd8_kw_xs_rank;
  KW_stat_cd4cd8=_KW_;
  KW_stat_cd4cd8_P=P_KW;
  if P_KW = . then flag_kw_cd4cd8_p05=.;
  else if P_KW < 0.05 then flag_kw_cd4cd8_p05=1;
  else flag_kw_cd4cd8_p05=0;
  keep gene_id KW_stat_cd4cd8 KW_stat_cd4cd8_P flag_kw_cd4cd8_p05;
run;

data cd4cd19_kw_2;
  set cd4cd19_kw_xs_rank;
  KW_stat_cd4cd19=_KW_;
  KW_stat_cd4cd19_P=P_KW;
  if P_KW = . then flag_kw_cd4cd19_p05=.;
  else if P_KW < 0.05 then flag_kw_cd4cd19_p05=1;
  else flag_kw_cd4cd19_p05=0;
  keep gene_id KW_stat_cd4cd19 KW_stat_cd4cd19_P flag_kw_cd4cd19_p05;
run;

data cd8cd19_kw_2;
  set cd8cd19_kw_xs_rank;
  KW_stat_cd8cd19=_KW_;
  KW_stat_cd8cd19_P=P_KW;
  if P_KW = . then flag_kw_cd8cd19_p05=.;
  else if P_KW < 0.05 then flag_kw_cd8cd19_p05=1;
  else flag_kw_cd8cd19_p05=0;
  keep gene_id KW_stat_cd8cd19 KW_stat_cd8cd19_P flag_kw_cd8cd19_p05;
run;

proc sort data=cd4cd8cd19_kw_2;
  by gene_id;
proc sort data=cd4cd8_kw_2;
  by gene_id;
proc sort data=cd4cd19_kw_2;
  by gene_id;
proc sort data=cd8cd19_kw_2;
  by gene_id;
run;


data kw_test_all;
  merge cd4cd8cd19_kw_2 cd4cd8_kw_2 cd4cd19_kw_2 cd8cd19_kw_2;
  by gene_id;
run;

data event.t1d_case_only_kw_test;
  set kw_test_all;
run;


proc freq data=kw_test_all noprint;
  tables flag_kw_all_p05*flag_kw_cd4cd8_p05*flag_kw_cd4cd19_p05*flag_kw_cd8cd19_p05 / out=kw_count;
proc print data=kw_count;
run;


/*
                            flag_kw_    flag_kw_
  flag_kw_     flag_kw_     cd4cd19_    cd8cd19_
   all_p05    cd4cd8_p05       p05         p05      COUNT

      0            0            0           0        5485
      0            0            0           1          20
      0            0            1           0          23
      0            0            1           1           2
      0            1            0           0           4
      1            0            0           1          10
      1            0            1           0           7
      1            0            1           1         188
      1            1            0           0           1
      1            1            0           1           5
      1            1            1           0          11
      1            1            1           1          54

total of 325 genes with MEI differences between cell types and/or subjects
of 65 that are really-really different

*/


