ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Export fusion counts for making heatmaps. I want to also include the number of times the exon is the MEE */


data genes_mei;
  set event.t1d_case_only_kw_test_3lvl_reduc;
  if flag_kw_cd4cd8_p05=1 or flag_kw_cd4cd19_p05=1 or flag_kw_cd4cd19_p05=1;
  sum_mei_diff=sum(flag_kw_cd4cd8_p05,flag_kw_cd4cd19_p05,flag_kw_cd4cd19_p05);
  keep gene_id sum_mei_diff;
run;

data genes_mee;
  set event.t1d_case_only_kw_exons_3level_v2;
  if flag_kw_cd4cd8_p05=1 or flag_kw_cd4cd19_p05=1 or flag_kw_cd4cd19_p05=1;
  sum_mee_diff=sum(flag_kw_cd4cd8_p05,flag_kw_cd4cd19_p05,flag_kw_cd4cd19_p05);
  keep gene_id sum_mee_diff;
run;

/* I am going to split these lists into: MEI-only, MEE-only, MEI and MEE */

proc sort data=genes_mei;
  by gene_id;
proc sort data=genes_mee;
  by gene_id;
run;

data genes_mei_mee;
  merge genes_mei (in=in1) genes_mee (in=in2);
  by gene_id;
  if in1 then flag_diff_mei=1; else flag_diff_mei=0;
  if in2 then flag_diff_mee=1; else flag_diff_mee=0;
run;

proc freq data=genes_mei_mee;
   tables flaG_diff_mei*flag_diff_mee;
run;

/*
 flag_diff_mei     flag_diff_mee

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    654 |    654
          |   0.00 |  94.92 |  94.92
          |   0.00 | 100.00 |
          |   0.00 |  95.75 |
 ---------+--------+--------+
        1 |      6 |     29 |     35
          |   0.87 |   4.21 |   5.08
          |  17.14 |  82.86 |
          | 100.00 |   4.25 |
 ---------+--------+--------+
 Total           6      683      689
              0.87    99.13   100.00


I will output the data for all 689 genes, but I am only going to plot the 29 genes first to examine differences
as these are the genes that are apparently most different between cell types
*/

data frag2gene;
  set hg19.hg19_exon_fragment_flagged;
  length specificity $4.;
  length fragment_name $30.;
  length fragment_id2 $30.;
  if flag_unique=1 then specificity="UNIQ";
  else if flag_constitutive=1 then specificity="CNST";
  else specificity="CMMN";
  fragment_id2=catx(":",chr,fragment_start,fragment_end);
  fragment_name=catt(fragment_id,"|",specificity," (",num_xscripts_per_fragment,")");
  keep fragment_start fragment_id2 fragment_name gene_id;
  rename fragment_id2=fragment_id;
run;

data gene_strand;
   set hg19.unique_info_fusions_si;
   where fusion_strand ne "" and gene_id ^? "|";
   keep fusion_strand gene_id;
run;

proc sort data=frag2gene nodup;
  by gene_id;
proc sort data=gene_strand nodup;
  by gene_id;
proc sort data=genes_mei_mee;
  by gene_id;
run;

data frag2keep;
   merge genes_mei_mee (in=in1) frag2gene (in=in2) gene_strand;
  by gene_id;
  if in1 and in2;
  drop flag_diff_mei flag_diff_mee;
run;


data frag_counts;
  set eventloc.t1d_case_only_fragment_counts;
  log_apn=log(apn+1);
  keep name cell_type fragment_id log_apn;
run;

proc sort data=frag2keep nodup;
  by fragment_id gene_id;
proc sort data=frag_counts;
  by fragment_id;
run;

data frag_counts2keep;
  merge frag2keep (in=in1) frag_counts (in=in2);
  by fragment_id;
  if in1 and in2;
run;

data design;
  set con.design_by_subject_new;
  if subject_id in ("M075","M048","M080") then delete;
  keep name subject_id cell_Type;
run;

proc sort data=frag_counts2keep;
  by name;
proc sort data=design;
  by name;
run;

data frag2gene_counts_w_key;
  merge design (in=in1) frag_counts2keep (in=in2);
  by name;
  if in1 and in2;
run;

proc sort data=frag2gene_counts_w_key;
   by gene_id fusion_strand subject_id fragment_start fragment_name cell_type;
run;

proc transpose data=frag2gene_counts_w_key out=counts_sbys;
   by  gene_id fusion_strand subject_id fragment_start fragment_name;
   var log_apn;
   id cell_type;
run;

data counts;
  set counts_sbys;
  drop _NAME_;
  rename cd4=logapn_cd4 cd8=logapn_cd8 cd19=logapn_cd19;
run;

proc sort data=counts;
  by gene_id fusion_strand fragment_start fragment_name;
proc means data=counts noprint;
  by gene_id fusion_strand fragment_start fragment_name;
  var logapn_CD19 logapn_CD4 logapn_CD8;
  output out=mean_counts mean=;
run;

data mean_counts2;
  set mean_counts;
  drop _TYPE_ _FREQ_;
run;

/* Count times where a fragment is the MEF */


proc sort data=frag2gene_counts_w_key;
   by gene_id cell_type subject_id descending log_apn;
run;

data flag_mee;
  set frag2gene_counts_w_key;
  by gene_id cell_type subject_id;
  if first.subject_id then flag_mee=1;
  else flag_mee=0;
run;

proc sort data=flag_mee;
  by gene_id fragment_name cell_type ;
proc means data=flag_mee noprint;
  by gene_id fragment_name cell_type ;
   var flag_mee;
  output out=sum_flag_mee sum=;
run;

proc transpose data=sum_flag_mee out=mee_sbys;
  by gene_id fragment_name;
  id cell_type;
  var flag_mee;
run;

data mee_sbys2;
  set mee_sbys;
  drop _NAME_;
  rename CD19=flag_mee_cd19 CD4=flag_mee_cd4 CD8=flag_mee_cd8;
run;


proc sort data=mee_sbys2;
   by fragment_name;
proc sort data=mean_counts2;
   by fragment_name;
run;

data export_counts_mee;
  merge mean_counts2 (in=in1) mee_sbys2 (in=in2);
   by fragment_name;
  if in1 and in2 then output;
  else if in1 then do;
     flag_mee_cd19=0;
     flag_mee_cd4=0;
     flag_mee_cd8=0;
     output; end;
run;

proc sort data=export_counts_mee;
  by gene_id fragment_start;
run;

data export_counts_mee2;
  retain gene_id fusion_strand fragment_name logapn_cd19 logapn_cd4 logapn_cd8 flag_mee_cd19 flag_mee_cd4 flag_mee_cd8;
  set export_counts_mee;
  drop fragment_start ;
  rename fusion_strand=strand;
run;


proc export data=export_counts_mee2 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_genes_mei_and_mee_fragments.csv"
   dbms=csv replace;
run;

*output gene list;

data genes1 genes2 genes3;
   set genes_mei_mee;
   if  flaG_diff_mei=1 and flag_diff_mee=1 then output genes1;
   else if flaG_diff_mei=1 then output genes2;
   else if  flaG_diff_mee=1 then output genes3;
   keep gene_id;
run;

proc sort data=genes1 nodup;
  by gene_id;
proc sort data=genes2 nodup;
  by gene_id;
proc sort data=genes3 nodup;
  by gene_id;
run;

proc export data=genes1 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_gene_list_3level_rank_mei_and_mee_reduced.csv"
   dbms=csv replace; putnames=no;
run;

proc export data=genes2 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_gene_list_3level_rank_mei_only_reduced.csv"
   dbms=csv replace; putnames=no;
run;

proc export data=genes3 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_gene_list_3level_rank_mee_only_reduced.csv"
   dbms=csv replace; putnames=no;
run;


