ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";

/* Take list of genes with DD exons and genes that are QDS and merge with MISO genes
   (1) DD genes vs QDS vs MISO -- what is annotated?
   (2) Of the genes that are annotated in both, what genes are expressed?
   (3) Of the expressed genes, what genes are sig with MISO? BF>10, BF>5

*/

/* (2) Of the genes that are annotated in both, what genes were analyzable with MISO? (ie, there was detection) */


data miso_se;
   set event.miso_events_gene_in_refseq;
   keep ens_gene_id event_name;
run;

data miso_analyzed;
   set event.miso_all_results_nsc_v_old_trim;
   where NSC_OLD_diff ne .;
   keep event_name;
run;

proc sort data=miso_se;
  by event_name;
proc sort data=miso_analyzed;
  by event_name;
run;

data miso_analyzed2ens;
  merge miso_se (in=in1) miso_analyzed (in=in2);
  by event_name;
  if in1 and in2;
run; *2530 events analyzed;

data miso_genes;
  set miso_analyzed2ens;
  keep ens_gene_id;
run;

proc sort data=miso_genes nodup;
  by ens_gene_id;
run;

data ens2refseq;
   set event.ensembl2refseq_gene_id;
   keep ens_gene_id gene_id;
run;

proc sort data=miso_genes nodup;
   by ens_gene_id;
proc sort data=ens2refseq nodup;
   by ens_gene_id gene_id;
run;

data miso_se2refseq_dtct;
   merge miso_genes (in=in1) ens2refseq (in=in2);
   by ens_gene_id;
   if in1 and in2 then flag_has_refseq=1;
   else if in1 then flag_has_refseq=0;
   if in1 then output;
run;  *1667 RefSeq genes;


*Genes that can be assessed for alternative splicing;

data genes2keep;
  set event.miso_refseq_exonskip_compare_qds;
  keep gene_id;
run;


proc sort data=genes2keep nodup;
   by gene_id;
proc sort data=miso_se2refseq_dtct;
   by gene_id;
run;

data miso2es_dtct;
  merge genes2keep (in=in1) miso_se2refseq_dtct (in=in2);
   by gene_id;
  if in2 then flag_has_miso_se_dtct=1; else flag_has_miso_se_dtct=0;
  if in1 then output;
run;

*MISO genes;
data genes_to_keep;
   set event.miso_refseq_exonskip_compare_qds;
   keep gene_id flag_has_multiple_exons flag_has_miso_se;
run; *2014 genes;

proc sort data=genes_to_keep nodup;
   by gene_id;
proc sort data=miso2es_dtct nodup;
   by gene_id;
run;

data miso2es_dtct2;
  merge miso2es_dtct (in=in1) genes_to_keep (in=in2);
  by gene_id;
  if in1 and in2;
  if gene_id="" then delete;
run;

proc freq data=miso2es_dtct2 noprint;
  tables flag_has_miso_se*flag_has_miso_se_dtct*flag_has_refseq*flag_has_multiple_exons / out=es_count;
run;

proc print data=es_count;
run;

/*
              flag_has_                 flag_has_
 flag_has_     miso_se_    flag_has_    multiple_
  miso_se        dtct        refseq       exons      COUNT    PERCENT

     1            0            .            1          94         .
     1            1            1            1         654       100


Overlap is 654
94 genes that are only in the DD/DQS set
0 that are only MISO

*/
/* Make permenant */

data event.miso_refseq_exnskp_cmpr_dtct_qds;
   set miso2es_dtct2;
run;

