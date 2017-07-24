ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Check that the set of genes with detected MISO SE are the same set detected  in Event analysis
  For this, I only want to look at genes that both a MISO event and are expressed */

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
  set event.miso_refseq_exonskip_compare;
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
   set event.miso_refseq_exonskip_Compare;
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

      1            0            .            0            3      .
      1            0            .            1          282      .
      1            1            1            0           14     1.2400
      1            1            1            1         1115    98.7600


MISO ONLY: 0
Overlap is 1115 (+14 if including single-exon genes)
RefSeq only: 13178 (from previous counts)
, or 282 genes with MISO, but none analyzed (+3 if including single-exon genes)

*/
/* Make permenant */

data event.miso_refseq_exonskip_cmpr_dtct;
   set miso2es_dtct2;
run;

 

