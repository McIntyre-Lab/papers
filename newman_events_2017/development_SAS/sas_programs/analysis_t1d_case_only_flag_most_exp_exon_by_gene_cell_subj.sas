/* For genes that had transcripts quantified with RSEM, flag the most expressed exon (MEE) */

ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Get genes with RSEM estimates */

data xscripts;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
  keep transcript_id;
run;

data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  keep gene_id transcript_id;
run;

proc sort data=xscripts nodup;
   by transcript_id;
proc sort data=xs2gene;
   by transcript_id;
run;

data genes2keep;
  merge xscripts (in=in1) xs2gene (in=in2);
  by transcript_id;
  if in1 and in2;
  keep gene_id;
run;

data fus2gene;
  set hg19.hg19_aceview_fusions_si_info;
  keep gene_id fusion_id;
run;

proc sort data=genes2keep nodup;
  by gene_id;
proc sort data=fus2gene nodup;
  by gene_id fusion_id;
run;

data fus2keep;
  merge genes2keep (in=in1) fus2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Get counts for fusions of interest */

data fus_counts;
   set con.fusion_q3_norm_data_all;
   keep fusion_id name log_q3_q3_apn;
run;

data design;
  set con.design_by_subject_new;
  if subject_id in ('M048','M075','M080') then delete;
*  if name =  '2009-PC-0221' then delete; *sample 75 cd8;
*  if name =  '2009-PC-0144' then delete; *sample 48 cd4;
*  if name =  '2009-PC-0236' then delete; *sample 80;
*  if name =  '2009-PC-0237' then delete; *sample 80;
*  if name =  '2009-PC-0235' then delete; *sample 80;
  keep name subject_id cell_type;
run;

proc sort data=fus_counts;
  by name;
proc sort data=design;
  by name;
run;

data counts_w_key;
  merge design (in=in1) fus_counts (in=in2);
  by name;
  if in1 and in2;
run;

proc sort data=counts_w_key;
  by fusion_id;
proc sort data=fus2keep;
  by fusion_id;
run;

data counts_w_key_subset;
  merge fus2keep (in=in1) counts_w_key (in=in2);
  by fusion_id;
  if in1 and in2;
run;

/* For each gene, flag the MEE by subject and cell type.
   I may also need to do this on exon means too, as an "overall" check */

proc sort data=counts_w_key_subset;
   by gene_id cell_type subject_id descending log_q3_q3_apn fusion_id;
run;

data flag_mee_by_subj_cell;
  set counts_w_key_subset;
  by gene_id cell_type subject_id;
  if first.subject_id then do;
        flag_mee=1;
        exon_rank=1; end;
  else if last.subject_id then do;
        flag_mee=0;
        exon_rank=3;
        end;
  else do;
        flag_mee=0;
        exon_rank=2;
        end;
run;

/* Transpose MEE flag, exon rank and apn */

proc sort data=flag_mee_by_subj_cell;
  by gene_id fusion_id subject_id cell_type;
run;

proc transpose data=flag_mee_by_subj_cell out=mee_by_cell_sbys;
   by gene_id fusion_id subject_id;
   id cell_type;
   var flag_mee;
run;

proc transpose data=flag_mee_by_subj_cell out=apn_by_cell_sbys;
   by gene_id fusion_id subject_id;
   id cell_type;
   var log_q3_q3_apn;
run;

proc transpose data=flag_mee_by_subj_cell out=rank_by_cell_sbys;
   by gene_id fusion_id subject_id;
   id cell_type;
   var exon_rank;
run;


data mee_by_cell_sbys2;
  set mee_by_cell_sbys;
  drop _NAME_;
  rename CD19=flag_mee_cd19 CD4=flag_mee_cd4 CD8=flag_mee_cd8;
run;

data apn_by_cell_sbys2;
  set apn_by_cell_sbys;
  drop _NAME_;
  rename CD19=log_q3_apn_cd19 CD4=log_q3_apn_cd4 CD8=log_q3_apn_cd8;
run;

data rank_by_cell_sbys2;
  set rank_by_cell_sbys;
  drop _NAME_;
  rename CD19=rank_cd19 CD4=rank_cd4 CD8=rank_cd8;
run;

proc sort data=mee_by_cell_sbys2;
  by gene_id fusion_id subject_id;
proc sort data=apn_by_cell_sbys2;
  by gene_id fusion_id subject_id;
proc sort data=rank_by_cell_sbys2;
  by gene_id fusion_id subject_id;
run;

data eventloc.t1d_exons_mee_rank_by_cell_subj;
  merge apn_by_cell_sbys2 mee_by_cell_sbys2 rank_by_cell_sbys2;
  by gene_id fusion_id subject_id;
run;

