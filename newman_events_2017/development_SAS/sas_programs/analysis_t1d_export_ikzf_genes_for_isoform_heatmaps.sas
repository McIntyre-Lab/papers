/* Export counts for genes that have significantly different MEIs/MEEs, based on the 3-level rankings */

ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

*I want to make "complete" and "reduced" IKZF plots, so keeping this for flagging only;

data reduced_xs;
  set eventloc.xs_reduced_by_tpm_subj_cell_v1;
  where flag_tpm_reduced=1 and flag_tpm_deleted_lt04=0;
  keep transcript_id ;
run;

data xs_counts;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
run;

data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  where gene_id ? "IKZF";
  keep transcript_id gene_id;
run;

proc sort data=xs_counts;
  by transcript_id;
proc sort data=reduced_xs nodup;
  by transcript_id;
proc sort data=xs2gene;
  by transcript_id;
run;

data xs2gene_counts;
  merge xs2gene (in=in1) xs_counts (in=in2) reduced_xs (in=in3);
  by transcript_id;
  if in1 and in2;
  if in3 then flag_reduced2keep=1; else flag_reduced2keep=0;
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
   by gene_id subject_id transcript_id flag_reduced2keep cell_type;
proc transpose data=xs2gene_counts_w_key out=counts_sbys;
   by gene_id subject_id transcript_id flag_reduced2keep;
   var tpm;
   id cell_type;
run;

data counts;
  set counts_sbys;
  drop _NAME_;
  rename cd4=tpm_cd4 cd8=tpm_cd8 cd19=tpm_cd19;
run;

proc sort data=counts;
  by gene_id transcript_id flag_reduced2keep;
proc means data=counts noprint;
  by gene_id transcript_id flag_reduced2keep;
  var TPM_CD19 TPM_CD4 TPM_CD8;
  output out=mean_counts mean=;
run;

data mean_counts2;
  set mean_counts;
  drop _TYPE_ _FREQ_;
run;

/* Count times where a transcript is the MEI */

data mei_all;
   set eventloc.t1d_xs_mei_rank_by_cell_subj_v2;
   where transcript_id ? "IKZF";
   keep transcript_id subject_id flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
run;

data mei_subset;
   set eventloc.t1d_xs_mei_rank_cell_sbj_reduced;
   where transcript_id ? "IKZF";
   keep transcript_id subject_id flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
run;


proc sort data=mei_all;
   by transcript_id;
proc means data=mei_all noprint;
   by transcript_id;
   var flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
   output out=mei_count_by_xs_all sum=;
run;

proc sort data=mei_subset;
   by transcript_id;
proc means data=mei_subset noprint;
   by transcript_id;
   var flag_mei_cd19 flag_mei_cd4 flag_mei_cd8;
   output out=mei_count_by_xs_subset sum=;
run;


data mei_count_by_xs_all2;
   set mei_count_by_xs_all;
   drop _TYPE_ _FREQ_;
run;


data mei_count_by_xs_subset2;
   set mei_count_by_xs_subset;
   drop _TYPE_ _FREQ_;
run;

proc sort data=mei_count_by_xs_all2;
   by transcript_id;
proc sort data=mei_count_by_xs_subset2;
   by transcript_id;
proc sort data=mean_counts2;
   by transcript_id;
run;

data export_counts_subset;
  merge mean_counts2 (in=in1) mei_count_by_xs_subset2 (in=in2);
  by transcript_id;
  if flag_reduced2keep=0 then delete;
  if in1 and in2 then output;
  else if in1 then do;
     flag_mei_cd19=0;
     flag_mei_cd4=0;
     flag_mei_cd8=0;
     output; end;
  drop flag_reduced2keep;
run;


data export_counts_all;
  merge mean_counts2 (in=in1) mei_count_by_xs_all2 (in=in2);
  by transcript_id;
  if in1 and in2 then output;
  else if in1 then do;
     flag_mei_cd19=0;
     flag_mei_cd4=0;
     flag_mei_cd8=0;
     output; end;
  drop flag_reduced2keep;
run;


proc export data=export_counts_subset outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_IKZFs_reduced_nolow.csv"
   dbms=csv replace;
run;

proc export data=export_counts_all outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_IKZFs_reduced_all.csv"
   dbms=csv replace;
run;


/* Now all IKZF genes */
data xs_counts;
  set eventloc.hg19_rsem_all_xscripts;
run;

data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  where gene_id ? "IKZF";
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
   by gene_id subject_id transcript_id  cell_type;
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

proc sort data=counts;
  by gene_id transcript_id;
proc means data=counts noprint;
  by gene_id transcript_id ;
  var TPM_CD19 TPM_CD4 TPM_CD8;
  output out=mean_counts mean=;
run;

data mean_counts2;
  set mean_counts;
     flag_mei_cd19="n.d.";
     flag_mei_cd4="n.d.";
     flag_mei_cd8="n.d.";
  drop _TYPE_ _FREQ_;
run;

proc export data=mean_counts2 outfile="!MCLAB/event_analysis/analysis_output/t1d_cases_IKZFs_all.csv"
   dbms=csv replace;
run;

