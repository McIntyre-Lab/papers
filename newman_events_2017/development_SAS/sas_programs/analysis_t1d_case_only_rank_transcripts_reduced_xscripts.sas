/* Now for my set of reduced transcripts I want to re-run the rank test then output the list of genes with really big diffs */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

/* Rank transcripts by TPM and perform K-W test pairwise by gene (CD4 vs CD8, CD4 vs CD19, CD8 vs CD19) */

/* Try: exclude genes where there is one isoform that is the MEI in 80% of samples: these are likely cases were there are only a few individuals that are different but in general the MEI is consistent */

data mei_flags;
  set event.t1d_flag_mei_cell_reduced_sbys;
  keep gene_id transcript_id subject_id flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
run;

proc sort data=mei_flags;
  by gene_id transcript_id;
proc means data=mei_flags noprint;
  by gene_id transcript_id;
  var flag_mei_cd19 flag_mei_cd8 flag_mei_cd4;
  output out=perc_mei_by_cell mean=;
run;

* ID genes to drop;

data genes2drop;
  set perc_mei_by_cell;
  if flag_mei_cd19 ge 0.8
     and flag_mei_cd8 ge 0.8
     and flag_mei_cd4 ge 0.8
     then output;
  keep gene_id;
run;

proc sort data=genes2drop nodup;
    by gene_id;
run;


/* Select transcripts */
data reduced_xs;
  set eventloc.xs_reduced_by_tpm_subj_cell_v1;
*  where flag_tpm_reduced=1 and flag_tpm_deleted_lt1=0;
  where flag_tpm_reduced=1 and flag_tpm_deleted_lt04=0;
*  set eventloc.xs_reduced_by_tpm_subj_cell_v2;
*  where flag_tpm_reduced=1;
  keep gene_id transcript_id;
run;

data check;
  set reduced_xs;
  keep gene_id;
run;

proc sort data=check nodup;
  by gene_id;
run;


proc sort data=reduced_xs;
  by gene_id;
proc sort data=genes2drop;
  by gene_id;
run;

data reduced_xs2;
  merge reduced_xs (in=in1) genes2drop (in=in2);
  by gene_id;
  if in2 then delete;
  keep transcript_id;
run;


data xs_counts;
  *set eventloc.xs_reduced_by_tpm_subj_cell_v1;
  keep subject_id library cell_type transcript_id gene_id tpm;
run;

*proc sort data=reduced_xs2 nodup;
proc sort data=reduced_xs nodup;
  by transcript_id;
run;

proc sort data=xs_counts;
  by transcript_id;
run;

data xs_counts_reduced nocount;
  merge reduced_xs (in=in1) xs_counts (in=in2);
  by transcript_id;
  if in1 and in2 then output xs_counts_reduced;
  else if in1 then output nocount;
run;

data design;
  set con.design_by_subject_new;
  if subject_id in ("M075","M048","M080") then delete;
  keep library subject_id cell_Type;
run;

proc sort data=xs_counts_reduced;
  by library;
proc sort data=design;
  by library;
run;

data xs2gene_counts_w_key;
  merge design (in=in1) xs_counts_reduced (in=in2);
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

/*
data event.t1d_case_only_kw_test_reduced;
  set kw_test_all;
run;
*/

proc freq data=kw_test_all noprint;
  tables flag_kw_all_p05*flag_kw_cd4cd8_p05*flag_kw_cd4cd19_p05*flag_kw_cd8cd19_p05 / out=kw_count;
proc print data=kw_count;
run;


/*
                           flag_kw_    flag_kw_
 flag_kw_     flag_kw_     cd4cd19_    cd8cd19_
  all_p05    cd4cd8_p05       p05         p05      COUNT

     0            0            0           0        5755
     0            0            0           1           2
     0            0            1           0           3
     1            0            1           1          15
     1            1            1           0           3
     1            1            1           1           9

total of 32 genes that have different MEIs after removing low-expressed transcripts
and 9 that are really-really different.

                             flag_kw_    flag_kw_
   flag_kw_     flag_kw_     cd4cd19_    cd8cd19_
    all_p05    cd4cd8_p05       p05         p05      COUNT    PERCENT

       0            0            0           0        5691    98.0700
       0            0            0           1           4     0.0689
       0            0            1           0           6     0.1034
       1            0            0           1           1     0.0172
       1            0            1           0           7     0.1206
       1            0            1           1          54     0.9306
       1            1            0           0           1     0.0172
       1            1            0           1           3     0.0517
       1            1            1           0           4     0.0689
       1            1            1           1          32     0.5514

*/
