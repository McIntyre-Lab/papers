libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

/* For the autoimmune genes that were differentially spliced in the case-only analysis of events
   I want to make a heatmap of transcripts per cell type visualize possible shifts in isoform expression

   I am only going to include the genes used in Event Analysis, and averaging RSEM estimates over individuals */

/* Get list of autoimmune genes */

data immunogenes;
  set con.immunobase_gene_flags;
  where flag_immuno_gene=1;
  keep gene_id;
run; *1854 possible genes;

data ds_genes;
  set con.flag_genes_dd_ds_de;
  where sum(flag_cd4_gene_on,flag_cd8_gene_on,flag_cd19_gene_on) > 1
    and flag_gene_ds=1;
  keep gene_id;
run;

proc sort data=immunogenes nodup;
  by gene_id;
proc sort data=ds_genes nodup;
  by gene_id;
run;

data genes2keep;
  merge ds_genes (in=in1) immunogenes (in=in2);
  by gene_id;
  if in1 and in2;
run;


/* Subset transcripts from reduced T1D set */

data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  keep transcript_id gene_id;
run;

proc sort data=genes2keep nodup;
   by gene_id;
proc sort data=xs2gene;
   by gene_id transcript_id;
run;

data xs2keep;
  merge genes2keep (in=in1) xs2gene (in=in2);
  by gene_id;
  if in1 and in2;
run; *22385 possible transcripts total;

/* Merge counts and calculate the cell-type average */

data xs_counts;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
run;

proc sort data=xs2keep;
  by transcript_id;
proc sort data=xs_counts;
  by transcript_id;
run;

data xs_counts_immuno;
  merge xs2keep (in=in1) xs_counts (in=in2);
  by transcript_id;
  if in1 and in2;
run; *790152 obs/246 samples = 3212 transcripts ;

data design;
  set con.design_by_subject_new;
  if name =  '2009-PC-0221' then delete; *sample 75 cd8;
  if name =  '2009-PC-0144' then delete; *sample 48 cd4;
  if name =  '2009-PC-0236' then delete; *sample 80;
  if name =  '2009-PC-0237' then delete; *sample 80;
  if name =  '2009-PC-0235' then delete; *sample 80; 
  keep library name cell_type;
run;

proc sort data=xs_counts_immuno;
  by library;
proc sort data=design;
  by library;
run;

data xs_counts_w_cell;
  merge design (in=in1) xs_counts_immuno (in=in2);
  by library;
  if in1 and in2;
run;

proc sort data=xs_counts_w_cell;
  by gene_id transcript_id cell_type;
proc means data=xs_counts_w_cell noprint;
  by gene_id transcript_id cell_type;
  var tpm;
  output out=mean_tpm_by_xs_cell mean=;
run;

data mean_tpm_by_xs_cell2;
  set mean_tpm_by_xs_cell;
  log_tpm=log(tpm+1);
run;

proc sort data=mean_tpm_by_xs_cell2;
   by gene_id transcript_id cell_type;
proc transpose data=mean_tpm_by_xs_cell2 out=log_tpm_by_xs_cell_sbys;
   by gene_id transcript_id;
   id cell_type;
   var log_tpm;
run;

data t1d_genes;
  set con.immunogene_flags;
  where flag_diabetes_gene=1;
  keep gene_id;
run;

proc sort data=t1d_genes nodup;
  by gene_id;
proc sort data=log_tpm_by_xs_cell_sbys;
  by gene_id transcript_id;
run;

data log_tpm_by_xs_cell_sbys2;
  merge t1d_genes (in=in1) log_tpm_by_xs_cell_sbys (in=in2);
  by gene_id;
  if in1 and in2;
run; *876 transcripts;

proc export data=log_tpm_by_xs_cell_sbys2 outfile="!MCLAB/event_analysis/analysis_output/mean_log_tpm_by_xscript_celltype_t1d_only.csv"
     dbms=csv replace;
run;


data ctla4 il7r ikzf1 ikzf3 ;
   set log_tpm_by_xs_cell_sbys2;
   if transcript_id="CTLA4.cAug10"
   or transcript_id="IKZF3.cAug10"
   or transcript_id="IKZF3.hAug10"
   or transcript_id="IKZF3.jAug10"
   or transcript_id="IKZF3.pAug10"
   then delete;
   if gene_id="CTLA4" then output ctla4;
   if gene_id="IL7R" then output  il7r;
   if gene_id="IKZF1" then output ikzf1;
   if gene_id="IKZF3" then output ikzf3;
run;

proc export data=ctla4
     outfile="!MCLAB/event_analysis/analysis_output/mean_log_tpm_by_xscript_celltype_t1d_only_CTLA4.csv"
     dbms=csv replace;
run;
proc export data=il7r
     outfile="!MCLAB/event_analysis/analysis_output/mean_log_tpm_by_xscript_celltype_t1d_only_IL7R.csv"
     dbms=csv replace;
run;
proc export data=ikzf1
     outfile="!MCLAB/event_analysis/analysis_output/mean_log_tpm_by_xscript_celltype_t1d_only_IKZF1.csv"
     dbms=csv replace;
run;
proc export data=ikzf3
     outfile="!MCLAB/event_analysis/analysis_output/mean_log_tpm_by_xscript_celltype_t1d_only_IKZF3.csv"
     dbms=csv replace;
run;




