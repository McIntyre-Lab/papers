ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname con '!PATCON/sas_data';

/* Reducing the number of transcripts per gene to only those that are appreciably expressed (TPM>0.1)

If gene has fewer than 3 transcripts, then eliminate none
Else:
   count the number of transcripts with TPM > 1
   if num_TPM1_xs le (mum_xs_per_gene - 2), and TPM ge 0.4 then put transcript into low_exp
   if num_TPM1_xs le (mum_xs_per_gene - 2), and TPM ge 1 then put transcript into tpm1_exp;

Be sure to do this by cell type too!!!

Then I want to see where each transcript falls: ideally, I only want to keep the "good ones" (i.e. TPM>1)

*/


/* Calc transcripts per gene */

data xs2gene;
  set event.hg19_aceview_xs2gene_fasta_index;
  keep transcript_id gene_id;
run;

data xs_list;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
  keep transcript_id;
run;

proc sort data=xs2gene nodup;
  by transcript_id gene_id ;
proc sort data=xs_list nodup;
  by transcript_id;
run;

data xs2gene2;
  merge xs2gene (in=in1) xs_list (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs2gene2;
  by gene_id transcript_id;
proc freq data=xs2gene2 noprint;
  tables gene_id / out=xs_per_gene;
run;

data xs_per_gene2;
  set xs_per_gene;
  keep gene_id count;
  rename count=xscripts_per_gene;
run;

/* For each sample and each gene, calculate the number of transcripts with TPM > 1 */

data flag_tpm1;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
  if tpm < 1 then flag_tpm_lt1=1;
  else flag_tpm_lt1=0;
run;

proc sort data=flag_tpm1;
  by transcript_id;
proc sort data=xs2gene2;
  by transcript_id;
run;

data flag_tpm1_w_gene;
  merge xs2gene2 (in=in1) flag_tpm1 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=flag_tpm1_w_gene;
  by gene_id library;
proc means data=flag_tpm1_W_gene noprint;
  by gene_id library;
  var flag_tpm_lt1;
  output out=tpm1_xs_per_gene sum(flag_tpm_lt1)=xs_per_gene_tpm_lt1;
run;

proc sort data=tpm1_xs_per_gene;
  by gene_id;
proc sort data=xs_per_gene2;
  by gene_id;
run;

data xs_tpm1_per_gene;
  merge xs_per_gene2 (in=in1) tpm1_xs_per_gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Merge this back onto counts data and start eliminating transcripts */

proc sort data=xs_tpm1_per_gene;
   by gene_id library;
proc sort data=flag_tpm1_w_gene;
   by gene_id library;
run;

data xs_tpm1_per_gene_2;
  merge xs_tpm1_per_gene (in=in1) flag_tpm1_w_gene (in=in2);
  by gene_id library;
  if in1 and in2;
run;

*all in cell type info;
data design;
  set con.design_by_subject_new;
  if subject_id in ("M075","M048","M080") then delete;
  keep subject_id library cell_type;
run;

proc sort data=design;
  by library;
proc sort data=xs_tpm1_per_gene_2;
  by library;
run;

data xs_tpm1_per_gene_3;
  merge design (in=in1) xs_tpm1_per_gene_2 (in=in2);
  by library;
  if in1 and in2;
run;



/* Eliminate transcripts */

data tpm_reduced tpm_deleted_04 tpm_deleted_10;
  set xs_tpm1_per_gene_3;
  if xscripts_per_gene le 3 then output tpm_reduced;
  else do;
      if (xs_per_gene_tpm_lt1 le (xscripts_per_gene - 2)) and tpm < 0.4 then output tpm_deleted_04;
      else if (xs_per_gene_tpm_lt1 le (xscripts_per_gene - 2)) and tpm < 1 then output tpm_deleted_10;
      else output tpm_reduced;
      end;
run;

data tpm_deleted_04_2;
   set tpm_deleted_04;
   keep cell_Type transcript_id;
run;

data tpm_deleted_10_2;
   set tpm_deleted_10;
   keep cell_Type transcript_id;
run;

data tpm_reduced_2;
   set tpm_reduced;
   keep cell_Type transcript_id;
run;

proc sort data=tpm_deleted_04_2 nodup;
  by cell_type transcript_id;
proc sort data=tpm_deleted_10_2 nodup;
  by cell_type transcript_id;
proc sort data=tpm_reduced_2 nodup;
  by cell_type transcript_id;
proc sort data=xs_tpm1_per_gene_3;
  by cell_type transcript_id;
run;

data xs_flag_xscript_bin;
  merge xs_tpm1_per_gene_3 (in=in1) tpm_reduced_2 (in=in2) tpm_deleted_04_2 (in=in3) tpm_deleted_10_2 (in=in4);
  by cell_type transcript_id;
  if in2 then flag_tpm_reduced=1; else flag_tpm_reduced=0;
  if in3 then flag_tpm_deleted_lt04=1; else flag_tpm_deleted_lt04=0;
  if in4 then flag_tpm_deleted_lt1=1; else flag_tpm_deleted_lt1=0;
  if in1 then output;
run;

data eventloc.xs_reduced_by_tpm_subj_cell_v1;
   set xs_flag_xscript_bin;
   drop _TYPE_ _FREQ_;
run;


/* Alternatively: calc mean TPM and eliminate

For each cell type, calculate the mean TPM, then eliminate transcripts based on mean count rather than by subject 

*/

data xs_counts;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
run;

data design;
  set con.design_by_subject_new;
  if subject_id in ("M075","M048","M080") then delete;
  keep subject_id library cell_type;
run;

proc sort data=design;
  by library;
proc sort data=xs_counts;
  by library;
run;

data xs_counts_w_key;
  merge design (in=in1) xs_counts (in=in2);
  by library;
  if in1 and in2;
run;

proc sort data=xs_counts_w_key;
   by transcript_id cell_type;
proc means data=xs_counts_w_key noprint;
   by transcript_id cell_type;
   var tpm;
   output out=mean_tpm_by_xs_cell mean=mean_tpm;
run;

data mean_tpm_by_xs_cell2;
  set mean_tpm_by_xs_cell;
  if mean_tpm < 1 then flag_tpm_lt1=1; else flag_tpm_lt1=0;
run;

proc sort data=mean_tpm_by_xs_cell2;
  by transcript_id;
proc sort data=xs2gene2;
  by transcript_id;
run;

data xs_mean_sbys_w_gene;
  merge xs2gene2 (in=in1) mean_tpm_by_xs_cell2 (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs_mean_sbys_w_gene; 
   by gene_id cell_type;
proc means data=xs_mean_sbys_w_gene noprint;
  by gene_id  cell_type;
  var flag_tpm_lt1;
  output out=num_xs_lt1
      sum(flag_tpm_lt1)=num_xs_tpm_lt1 ;
run;

proc sort data=num_xs_lt1;
  by gene_id;
proc sort data=xs_per_gene2;
  by gene_id;

data xs_per_gene3;
  merge xs_per_gene2 (in=in1) num_xs_lt1 (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=xs_per_gene3;
  by gene_id cell_type;
proc sort data=xs_mean_sbys_w_gene;
  by gene_id cell_type;
run;

data xs_mean_sbys_w_cnts;
  merge xs_mean_sbys_w_gene (in=in1) xs_per_gene3 (in=in2);
  by gene_id cell_type;
  if in1 and in2 ;
  drop _TYPE_ _FREQ_;
run;

data tpm_reduced tpm_delete_cd4 tpm_delete_cd8 tpm_delete_cd19;
   set xs_mean_sbys_w_cnts;
   if xscripts_per_gene le 3 then output tpm_reduced;
   else do;
      if (num_xs_tpm_lt1 le (xscripts_per_gene - 2)) and mean_tpm < 0.4 then do;
            if cell_type="CD4" then output tpm_delete_cd4;
            else if cell_type="CD8" then output tpm_delete_cd8;
            else if cell_type="CD19" then  output tpm_delete_cd19;
            end;
      else output tpm_reduced;
      end;
run;

data tpm_reduced_2;
  set tpm_reduced;
  keep transcript_id;
run;

data tpm_delete_cd4_2;
  set tpm_delete_cd4;
  keep transcript_id;
run;

data tpm_delete_cd8_2;
  set tpm_delete_cd8;
  keep transcript_id;
run;

data tpm_delete_cd19_2;
  set tpm_delete_cd19;
  keep transcript_id;
run;

proc sort data=xs_tpm1_per_gene_3;
  by transcript_id;
proc sort data=tpm_reduced_2 nodup;
  by transcript_id;
proc sort data=tpm_delete_cd4_2 nodup;
  by transcript_id;
proc sort data=tpm_delete_cd8_2  nodup;
  by transcript_id;
proc sort data=tpm_delete_cd19_2  nodup;
  by transcript_id;
run;


data tpm_flag_counts_to_keep;
  merge xs_tpm1_per_gene_3 (in=in1) tpm_reduced_2 (in=in2) tpm_delete_cd4_2 (in=in3)
           tpm_delete_cd8_2 (in=in4) tpm_delete_cd19_2 (in=in5) ;
   by transcript_id;
  if in2 then flag_tpm_reduced=1; else flag_tpm_reduced=0;
  if in3 then flag_tpm_deleted_cd4=1; else flag_tpm_deleted_cd4=0;
  if in4 then flag_tpm_deleted_cd8=1; else flag_tpm_deleted_cd8=0;
  if in5 then flag_tpm_deleted_cd19=1; else flag_tpm_deleted_cd19=0;
  if in1 then output;
    drop _TYPE_ _FREQ_ xs_per_gene_tpm_lt1 flag_tpm_lt1;
run;

data eventloc.xs_reduced_by_tpm_subj_cell_v2;
   set tpm_flag_counts_to_keep;
run;


