/* "Three level" rank tests for isoforms and exons */

ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';


/* Un-transpose data so that ranks are stacked and not SbyS */

data exon_rank;
  set eventloc.t1d_exons_mee_rank_cell_subj_v2;
  keep gene_id fusion_id subject_id rank_cd19 rank_cd8 rank_Cd4;
run;

proc transpose data=exon_rank out=exon_rank_stack;
   by gene_id fusion_id subject_id;
   var rank_cd19 rank_cd8 rank_cd4;
run;

data exon_rank_stack2;
   length cell_type $4.;
   set exon_rank_stack;
   if _NAME_ = "rank_cd19" then cell_type="CD19";
   if _NAME_ = "rank_cd8" then cell_type="CD8";
   if _NAME_ = "rank_cd4" then cell_type="CD4";
   rename COL1=exon_rank;
   drop _NAME_;
run;

proc sort data=exon_rank_stack2;
  by gene_id cell_type exon_rank;
run;

proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  class cell_type;
  var exon_rank;
  output out=cd4cd8cd19_kw_exon_rank WILCOXON;
run;

proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  where cell_type ne "CD8";
  class cell_type;
  var exon_rank;
  output out=cd4cd19_kw_exon_rank WILCOXON;
run;


proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  where cell_type ne "CD4";
  class cell_type;
  var exon_rank;
  output out=cd8cd19_kw_exon_rank WILCOXON;
run;


proc npar1way data=exon_rank_stack2 noprint ;
  by gene_id;
  where cell_type ne "CD19";
  class cell_type;
  var exon_rank;
  output out=cd4cd8_kw_exon_rank WILCOXON;
run;


data cd4cd8cd19_kw_2;
  set cd4cd8cd19_kw_exon_rank;
  KW_stat_all=_KW_;
  KW_stat_all_P=P_KW;
  if P_KW = . then flag_kw_all_p05=.;
  else if P_KW < 0.05 then flag_kw_all_p05=1;
  else flag_kw_all_p05=0;
  keep gene_id KW_stat_all KW_stat_all_P flag_kw_all_p05;
run;

data cd4cd8_kw_2;
  set cd4cd8_kw_exon_rank;
  KW_stat_cd4cd8=_KW_;
  KW_stat_cd4cd8_P=P_KW;
  if P_KW = . then flag_kw_cd4cd8_p05=.;
  else if P_KW < 0.05 then flag_kw_cd4cd8_p05=1;
  else flag_kw_cd4cd8_p05=0;
  keep gene_id KW_stat_cd4cd8 KW_stat_cd4cd8_P flag_kw_cd4cd8_p05;
run;

data cd4cd19_kw_2;
  set cd4cd19_kw_exon_rank;
  KW_stat_cd4cd19=_KW_;
  KW_stat_cd4cd19_P=P_KW;
  if P_KW = . then flag_kw_cd4cd19_p05=.;
  else if P_KW < 0.05 then flag_kw_cd4cd19_p05=1;
  else flag_kw_cd4cd19_p05=0;
  keep gene_id KW_stat_cd4cd19 KW_stat_cd4cd19_P flag_kw_cd4cd19_p05;
run;

data cd8cd19_kw_2;
  set cd8cd19_kw_exon_rank;
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

data event.t1d_case_only_kw_exons_3level_v2;
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

     0            0            0           0        5040
     0            0            0           1          53
     0            0            1           0          51
     0            0            1           1           7
     0            1            0           0           4
     0            1            1           0           1
     1            0            0           1          34
     1            0            1           0          24
     1            0            1           1         391
     1            1            0           0           5
     1            1            0           1          18
     1            1            1           0          43
     1            1            1           1         139
*/




