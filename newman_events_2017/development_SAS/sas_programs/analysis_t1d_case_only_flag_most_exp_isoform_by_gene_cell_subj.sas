/* For genes that had transcripts quantified with RSEM, flag the most expressed exon (MEE) */

ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

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
   by  gene_id cell_type subject_id descending tpm transcript_id;
run;

/* For each gene, flag the MEI by subject and cell type. */

data flag_mei_by_subj_cell;
  set xs2gene_counts_w_key;
  by gene_id cell_type subject_id;
  if first.subject_id then do;
        flag_mei=1;
        xscript_rank=1; end;
  else if last.subject_id then do;
        flag_mei=0;
        xscript_rank=3;
        end;
  else do;
        flag_mei=0;
        xscript_rank=2;
        end;
run;

/* Transpose MEI flag, exon rank and apn */

proc sort data=flag_mei_by_subj_cell;
  by gene_id transcript_id subject_id cell_type;
run;

proc transpose data=flag_mei_by_subj_cell out=mei_by_cell_sbys;
   by gene_id transcript_id subject_id;
   id cell_type;
   var flag_mei;
run;

proc transpose data=flag_mei_by_subj_cell out=tpm_by_cell_sbys;
   by gene_id transcript_id subject_id;
   id cell_type;
   var tpm;
run;

proc transpose data=flag_mei_by_subj_cell out=rank_by_cell_sbys;
   by gene_id transcript_id subject_id;
   id cell_type;
   var xscript_rank;
run;


data mei_by_cell_sbys2;
  set mei_by_cell_sbys;
  drop _NAME_;
  rename CD19=flag_mei_cd19 CD4=flag_mei_cd4 CD8=flag_mei_cd8;
run;

data tpm_by_cell_sbys2;
  set tpm_by_cell_sbys;
  drop _NAME_;
  rename CD19=tpm_cd19 CD4=tpm_cd4 CD8=tpm_cd8;
run;

data rank_by_cell_sbys2;
  set rank_by_cell_sbys;
  drop _NAME_;
  rename CD19=rank_cd19 CD4=rank_cd4 CD8=rank_cd8;
run;

proc sort data=mei_by_cell_sbys2;
  by gene_id transcript_id subject_id;
proc sort data=tpm_by_cell_sbys2;
  by gene_id transcript_id subject_id;
proc sort data=rank_by_cell_sbys2;
  by gene_id transcript_id subject_id;
run;

data eventloc.t1d_xs_mei_rank_by_cell_subj;
  merge tpm_by_cell_sbys2 mei_by_cell_sbys2 rank_by_cell_sbys2;
  by gene_id transcript_id subject_id;
run;

