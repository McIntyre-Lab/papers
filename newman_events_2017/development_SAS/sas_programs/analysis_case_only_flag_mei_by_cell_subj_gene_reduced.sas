/* For each cell type, subject and gene, flag the most expressed transcript. Then we want to test the agreement between subjects within a subject, and the agreement between cell types within a subject */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

* count number of XS per gene;
data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  keep gene_id transcript_id;
run;

data xs2keep;
  set eventloc.xs_reduced_by_tpm_subj_cell_v1;
  where flag_tpm_reduced=1 and flag_tpm_deleted_lt04=0;
  keep transcript_id;
run;

proc sort data=xs2gene;
  by transcript_id;
proc sort data=xs2keep nodup;
  by transcript_id;
run;

data xs2gene2;
  merge xs2gene (in=in1) xs2keep (in=in2);
  by transcript_id;
  if in1 and in2;
run;


proc freq data=xs2gene2 noprint;
   tables gene_id / out=xs_count;
run;

data flag_gene_multi_xs;
  set xs_count;
*  if count > 1 then flag_multi_xs=1;
*  else flag_multi_xs=0;
*  keep gene_id flag_multi_xs;
  if count > 1 then output;
  keep gene_id;
run;


proc sort data=xs2gene2;
  by gene_id;
proc sort data=flag_gene_multi_xs;
   by gene_id;
run;

data xs2gene_flag_multi;
  merge xs2gene2 (in=in1) flag_gene_multi_xs (in=in2);
   by gene_id;
 if in1 and in2;
run;


data xs_counts;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
run;


data design;
  set con.design_by_subject_new;
  /* I am going to drop subjects if any data are missing */
*if name =  '2009-PC-0221' then delete; *sample 75 cd8;
*if name =  '2009-PC-0144' then delete; *sample 48 cd4;
*if name =  '2009-PC-0236' then delete; *sample 80;
*if name =  '2009-PC-0237' then delete; *sample 80;
*if name =  '2009-PC-0235' then delete; *sample 80; 
  if subject_id in ("M075","M048","M080") then delete;
  keep library subject_id Name cell_type;
run;

proc sort data=xs_counts;
  by library;
proc sort data=design;
  by library;
run;

data xs_counts_w_key;
  merge design (in=in1) xs_counts (in=in2);
  by library;
  if in1 and in2;
run;

proc sort data=xs_counts_w_key;
   by transcript_id;
proc sort data=xs2gene_flag_multi;
   by transcript_id;
run;

data xs2gene_w_key;
  merge xs2gene_flag_multi (in=in1) xs_counts_w_key (in=in2);
  by transcript_id;
  if in1 and in2;
run;


/* Now drop transcripts that are not in the reduced set */

proc sort data=xs2keep;
   by transcript_id;
proc sort data=xs2gene_w_key;
   by transcript_id;
run;

data flag_xs_low_exp;
  merge xs2keep (in=in1) xs2gene_w_key (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* Sort by cell type, gene, subject, and descending TPM */

proc sort data=flag_xs_low_exp;
   by cell_type gene_id subject_id descending TPM;
run;

data flag_mei;
   set flag_xs_low_exp;
   by cell_type gene_id subject_id;
   if first.subject_id then flag_mei=1;
   else flag_mei=0;
run;

proc sort data=flag_mei;
   by gene_id subject_id transcript_id cell_type;
proc transpose data=flag_mei out=flag_mei_cell_sbys;
   by gene_id subject_id transcript_id;
   id cell_type;
   var flag_mei;
run;
proc transpose data=flag_mei out=tpm_cell_sbys;
   by gene_id subject_id transcript_id;
   id cell_type;
   var tpm;
run;

data flag_mei_cell_sbys2;
   set flag_mei_cell_sbys;
   drop _NAME_;
   rename CD19=flag_mei_CD19 CD4=flag_mei_CD4 CD8=flag_mei_CD8;
run;

data tpm_cell_sbys2;
   set tpm_cell_sbys;
   drop _NAME_;
  rename CD19=TPM_CD19 CD4=TPM_CD4 CD8=TPM_CD8;
run;

proc sort data=flag_mei_cell_sbys2;
   by gene_id subject_id transcript_id;
proc sort data=tpm_cell_sbys2;
   by gene_id subject_id transcript_id;
run;


/* Make permenant */
     
data event.t1d_flag_mei_cell_reduced_sbys;
  merge flag_mei_cell_sbys2 (in=in1) tpm_cell_sbys2 (in=in2);
  by gene_id subject_id transcript_id;
  if in1 and in2;
run;

proc sort data=event.t1d_flag_mei_cell_reduced_sbys;;
   by gene_id subject_id;
run;

proc freq data=event.t1d_flag_mei_cell_reduced_sbys noprint;
    by gene_id;
    tables flag_mei_CD4*flag_mei_CD8 ;
    test agree;
    output out=cd4_cd8_agree AGREE;
run;


proc freq data=event.t1d_flag_mei_cell_reduced_sbys  noprint;
    by gene_id;
    tables flag_mei_CD4*flag_mei_CD19 ;
    test agree;
    output out=cd4_cd19_agree AGREE;
run;


proc freq data=event.t1d_flag_mei_cell_reduced_sbys  noprint;
    by gene_id ;
    tables flag_mei_CD8*flag_mei_CD19 ;
    test agree;
    output out=cd8_cd19_agree AGREE;
run;

data event.t1d_mei_by_cell_cd48_reduced;
   set cd4_cd8_agree;
run;
data event.t1d_mei_by_cell_cd419_reduced;
   set cd4_cd19_agree;
run;
data event.t1d_mei_by_cell_cd8d19_reduced;
   set cd8_cd19_agree;
run;

