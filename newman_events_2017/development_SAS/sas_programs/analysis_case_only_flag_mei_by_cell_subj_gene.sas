/* For each cell type, subject and gene, flag the most expressed transcript. Then we want to test the agreement between subjects within a subject, and the agreement between cell types within a subject */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
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


/* Now drop transcripts that are lowly-expressed */

data flag_xs_low_exp;
   set xs2gene_w_key;
   *if flag_multi_xs=1 and TPM < 1 then delete;
run;

/* Sort by cell type, gene, subject, and descending TPM */

proc sort data=xs2gene_w_key;
   by cell_type gene_id subject_id descending TPM;
run;

data flag_mei;
   set xs2gene_w_key;
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

/* Make permenant */

data event.t1d_flag_mei_cell_sbys;
  set flag_mei_cell_sbys;
  drop _NAME_;
run;



/* Before testing I want to have an idea as to the number of genes that have different MEI's between cell types */

data flag_mei_diff;
  set flag_mei_cell_sbys;
  *if CD4=0 and CD8=0 and CD19=0 then delete; *remove observations without MEI;
  if CD4=. or CD8=. or CD19=. then delete; *remove observations with missing;
  if CD4 ne CD8 then flag_CD4_CD8_MEI_diff=1;
  else flag_CD4_CD8_MEI_diff=0;
  if CD4 ne CD19 then flag_CD4_CD19_MEI_diff=1;
  else flag_CD4_CD19_MEI_diff=0;
  if CD8 ne CD19 then flag_CD8_CD19_MEI_diff=1;
  else flag_CD8_CD19_MEI_diff=0;
run;

proc sort data=flag_mei_diff;
  by gene_id;
proc means data=flag_mei_diff noprint;
   by gene_id;
  var flag_CD4_CD8_MEI_diff flag_CD4_CD19_MEI_diff flag_CD8_CD19_MEI_diff;
  output out=flag_mei_diff_gene max=;
run;

proc freq data=flag_mei_diff_gene;
  tables flag_CD4_CD8_MEI_diff flag_CD4_CD19_MEI_diff flag_CD8_CD19_MEI_diff;
run;



/* Check: per subject and gene, is there agreement with the MEI (pairwise cell types) */

/* across all samples:

    flag_CD4_CD8_                             Cumulative    Cumulative
         MEI_diff    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        3454       59.45          3454        59.45
                1        2356       40.55          5810       100.00


   flag_CD4_CD19_                             Cumulative    Cumulative
         MEI_diff    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        3667       63.12          3667        63.12
                1        2143       36.88          5810       100.00


   flag_CD8_CD19_                             Cumulative    Cumulative
         MEI_diff    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        3393       58.40          3393        58.40
                1        2417       41.60          5810       100.00

So I except ~40% of genes to have a different MEI

After dropping 

       flag_CD4_CD8_                             Cumulative    Cumulative
            MEI_diff    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0         627       21.02           627        21.02
                   1        2356       78.98          2983       100.00


      flag_CD4_CD19_                             Cumulative    Cumulative
            MEI_diff    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0         840       28.16           840        28.16
                   1        2143       71.84          2983       100.00


      flag_CD8_CD19_                             Cumulative    Cumulative
            MEI_diff    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0         566       18.97           566        18.97
                   1        2417       81.03          2983       100.00

 */

proc sort data=event.t1d_flag_mei_cell_sbys;;
   by gene_id subject_id;
run;

proc freq data=event.t1d_flag_mei_cell_sbys noprint;
    by gene_id;
    tables CD4*CD8 ;
    test agree;
    output out=cd4_cd8_agree AGREE;
run;


proc freq data=event.t1d_flag_mei_cell_sbys  noprint;
    by gene_id;
    tables CD4*CD19 ;
    test agree;
    output out=cd4_cd19_agree AGREE;
run;


proc freq data=event.t1d_flag_mei_cell_sbys  noprint;
    by gene_id ;
    tables CD8*CD19 ;
    test agree;
    output out=cd8_cd19_agree AGREE;
run;

