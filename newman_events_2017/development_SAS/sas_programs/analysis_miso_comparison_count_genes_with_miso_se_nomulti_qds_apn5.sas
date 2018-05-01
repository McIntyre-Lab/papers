ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";

/* Take list of genes with DD exons and genes that are QDS and merge with MISO genes
   (1) DD genes vs QDS vs MISO -- what is annotated?
   (2) Of the genes that are annotated in both, what genes are expressed?
   (3) Of the expressed genes, what genes are sig with MISO? BF>10, BF>5

*/

/* (1) What DD/QDS genes also have MISO annotations ? */

data miso_se;
   set event.miso_se2gene;
   keep ens_gene_id;
run;

data ens2refseq;
   set event.ensembl2refseq_gene_id;
   keep ens_gene_id gene_id;
run;

data genes2keep;
   set event.mm10_flag_gene_dd_ds_exons_apn5;
   if flag_gene_cell_specific ne 0 then delete;
   if flag_gene_monoexon=1 and flag_cell_by_fus_fdr05=. then delete;
   keep gene_id;
run;

proc sort data=ens2refseq nodup;
   by gene_id;
proc sort data=genes2keep nodup;
   by gene_id;
run;

data ens2refseq2 noens;
  merge ens2refseq (in=in1) genes2keep (in=in2);
  by gene_id;
  if in1 and in2 then output ens2refseq2;
  else if in2 then output noens;
run; *6328 of 20875 with ENSG match;
 * this is actually 6171 genes with at least one ENSG match, and 101 without a match;

proc sort data=miso_se nodup;
   by ens_gene_id;
proc sort data=ens2refseq2 nodup;
   by ens_gene_id gene_id;
run;

data miso_se2refseq no_miso;
   merge miso_se (in=in1) ens2refseq2 (in=in2);
   by ens_gene_id;
   if in1 and in2 then flag_has_refseq=1;
   else if in1 then flag_has_refseq=0;
   if in1 then output miso_se2refseq;
   else if in2 then output no_miso;
run;

data no_miso2;
  set no_miso;
  keep gene_id;
run;

proc sort data=no_miso2 nodup;
  by gene_id;
run; *5446 unique genes without MISO;


proc freq data=miso_se2refseq;
  tables flag_has_refseq;
run;

/*

                                            Cumulative    Cumulative
flag_has_refseq    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0        1788       70.50          1788        70.50
              1         748       29.50          2536       100.00

+ 101 genes tested without an ENSG match, and 5446 genes without a MISO match;

*/

/* Genes examined for DD or QDS */

*Genes with at least 2 exons;

data fus2gene;
   set mm10.mm10_refseq_fusion_si_info_v2;
   keep fusion_id primary_gene_id;
   rename primary_gene_id=gene_id;
run;

proc sort data=fus2gene nodup;
   by fusion_id gene_id;
proc freq data=fus2gene noprint;
   tables gene_id / out=fus_count_per_gene;
run;

data gene_2fus;
  set fus_count_per_gene;
  where count > 1;
  keep gene_id;
run;

proc sort data=miso_se2refseq nodup;
   by gene_id;
proc sort data=gene_2fus nodup;
   by gene_id;
run;

data miso2es;
  merge miso_se2refseq (in=in1) gene_2fus (in=in2);
   by gene_id;
  if in1 then flag_has_miso_se=1; else flag_has_miso_se=0;
  if in2 then flag_has_multiple_exons=1; else flag_has_multiple_exons=0;
  if in1 then output;
run;

proc freq data=miso2es;
  tables flag_has_miso_se*flag_has_refseq;
run;

proc freq data=miso2es;
  where flag_has_multiple_exons=1;
  tables flag_has_miso_se*flag_has_refseq;
run;

proc freq data=miso2es;
  where flag_has_refseq=1;
  tables flag_has_miso_se*flag_has_multiple_exons;
run;


/*
flag_has_miso_se
          flag_has_refseq

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       1 |   1788 |    748 |   2536
         |  70.50 |  29.50 | 100.00
         |  70.50 |  29.50 |
         | 100.00 | 100.00 |
---------+--------+--------+
Total        1788      748     2536
            70.50    29.50   100.00

 flag_has_miso_se
           flag_has_refseq

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        1 |    748 |    748
          | 100.00 | 100.00
          | 100.00 |
          | 100.00 |
 ---------+--------+
 Total         748      748
            100.00   100.00

   flag_has_miso_se
             flag_has_multiple_exons

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       1|  Total
   ---------+--------+
          1 |    748 |    748
            | 100.00 | 100.00
            | 100.00 |
            | 100.00 |
   ---------+--------+
   Total         748      748
              100.00   100.00



*/

data single_fusion;
   set miso2es;
   where flag_has_multiple_exons=0 and flag_has_refseq=1 and flag_has_miso_se=1;
run; *01 genes that have only a single exonic region, but have a MISO skipped exon;


/* make permenant */

data event.miso_refseq_exonskip_compare_qds;
   set miso2es;
run;




gedit analysis_miso_comparison_count_genes_with_miso_se_nomulti.sas
gedit analysis_miso_comparison_count_genes_with_analyzed_miso_se_jrbn2_nomulti.sas
gedit analysis_miso_comparison_count_genes_with_DE_miso_se_nomulti.sas

