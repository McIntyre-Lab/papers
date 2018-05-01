/* For each cell type, subject and gene, flag the most expressed transcript. Then we want to test the agreement between subjects within a subject, and the agreement between cell types within a subject */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
*libname event 'S:\McIntyre_Lab\event_analysis\sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

/* First for each gene, I want to remove transcripts with counts < 1, but only if there are multiple transcripts for a gene. If a gene has only one transcript, then I want to keep it

I also only want to do this for my reduced transcript list so first I need to filter transcripts */

* count number of XS per gene;
data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  keep gene_id transcript_id;
run;

data xs2keep;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
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
  if count > 1 then flag_mono_isoform_gene=0; else flag_mono_isoform_gene=1;
  keep gene_id flag_mono_isoform_gene;
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

/* Calc mean TPM */

proc sort data=xs2gene_w_key;
   by gene_id transcript_id flag_mono_isoform_gene cell_type;
proc means data=xs2gene_w_key noprint;
   by gene_id transcript_id flag_mono_isoform_gene cell_type;
   var tpm;
   output out=xs2gene_means mean=;
run;

/* Sort by cell type, gene, subject, and descending TPM */

proc sort data=xs2gene_means;
   by gene_id cell_type descending TPM;
run;

data flag_mei;
   set xs2gene_means;
   by gene_id cell_type;
   if first.cell_type then do;
       flag_mei=1;
       isoform_rank=1;
       end;
   else if last.cell_type then do;
       flag_mei=0;
       isoform_rank=3;
       end;
   else do;
       flag_mei=0;
       isoform_rank=2;
       end;
run;

proc sort data=flag_mei;
   by gene_id transcript_id flag_mono_isoform_gene cell_type;
proc transpose data=flag_mei out=flag_mei_cell_sbys;
   by gene_id transcript_id flag_mono_isoform_gene;
   id cell_type;
   var flag_mei;
run;
proc transpose data=flag_mei out=rank_cell_sbys;
   by gene_id transcript_id ;
   id cell_type;
   var isoform_rank;
run;
proc transpose data=flag_mei out=tpm_cell_sbys;
   by gene_id transcript_id;
   id cell_type;
   var tpm;
run;

data flag_mei_cell_sbys2;
   set flag_mei_cell_sbys;
   drop _NAME_;
   rename CD19=flag_mei_CD19 CD4=flag_mei_CD4 CD8=flag_mei_CD8;
run;

data rank_cell_sbys2;
   set rank_cell_sbys;
   drop _NAME_;
   rename CD19=isoform_rank_CD19 CD4=isoform_rank_CD4 CD8=isoform_rank_CD8;
run;

data tpm_cell_sbys2;
   set tpm_cell_sbys;
   drop _NAME_;
  rename CD19=TPM_CD19 CD4=TPM_CD4 CD8=TPM_CD8;
run;

proc sort data=flag_mei_cell_sbys2;
   by gene_id transcript_id;
proc sort data=rank_cell_sbys2;
   by gene_id transcript_id;
proc sort data=tpm_cell_sbys2;
   by gene_id transcript_id;
run;


/* Make permenant */

data event.t1d_flag_mei_cell_sbys_mean;
  merge flag_mei_cell_sbys2 (in=in1) tpm_cell_sbys2 (in=in2) rank_cell_sbys2 (in=in3);
  by gene_id transcript_id;
  if in1 and in2 and in3;
run;


/* Now count genes with MEI differences */

data cd4 cd8 cd19;
   set event.t1d_flag_mei_cell_sbys_mean;
   if flag_mei_cd19=1 then output cd19;
   if flag_mei_cd4=1 then output cd4;
   if flag_mei_cd8=1 then output cd8;
   keep gene_id transcript_id flag_mono_isoform_gene;
run;

data cd4_2;
  set cd4;
  rename transcript_id=cd4_mei_transcript_id;
run;

data cd8_2;
  set cd8;
  rename transcript_id=cd8_mei_transcript_id;
  drop flag_mono_isoform_gene;
run;

data cd19_2;
  set cd19;
  rename transcript_id=cd19_mei_transcript_id;
  drop flag_mono_isoform_gene;
run;

proc sort data=cd4_2;
  by gene_id;
proc sort data=cd19_2;
  by gene_id;
proc sort data=cd8_2;
  by gene_id;
run;

data flag_gene_mei_diff;
  merge cd4_2 cd8_2 cd19_2;
  by gene_id;
run;

data flag_gene_mei_diff2;
  set flag_gene_mei_diff;
  if cd4_mei_transcript_id=cd8_mei_transcript_id then flag_cd4cd8_mei_diff=0; else flag_cd4cd8_mei_diff=1;
  if cd4_mei_transcript_id=cd19_mei_transcript_id then flag_cd4cd19_mei_diff=0; else flag_cd4cd19_mei_diff=1;
  if cd8_mei_transcript_id=cd19_mei_transcript_id then flag_cd8cd19_mei_diff=0; else flag_cd8cd19_mei_diff=1;
run;

proc freq data=flag_gene_mei_diff2 noprint;
 tables flag_cd4cd8_mei_diff*flag_cd4cd19_mei_diff*flag_cd8cd19_mei_diff*flag_mono_isoform_gene / out=mei_diff_check;
proc print data=mei_diff_check;
run;
quit;

/*
                                                   flag_mono_
 flag_cd4cd8_    flag_cd4cd19_    flag_cd8cd19_     isoform_
   mei_diff         mei_diff         mei_diff         gene       COUNT

       0               0                0               0         2722
       0               0                0               1         2827
       0               1                1               0          169
       1               0                1               0           41
       1               1                0               0           44
       1               1                1               0            7

*/

/* Make permenant -- I want to look at this later in comparison with the MEIs */

data event.t1d_mei_comparison_by_gene;
  set flag_gene_mei_diff2;
run;
