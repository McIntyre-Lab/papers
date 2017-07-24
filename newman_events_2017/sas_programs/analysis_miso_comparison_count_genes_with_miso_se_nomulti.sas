/* MISO analysis:
1. Check that the set of genes with MISO SE are the same set with SE in Event analysis
2. Check that the set of genes with detected MISO SE are the same set with detected SE in Event analysis
3. Check that the set of genes with DE MISO SE are the same set with DE/DD SE in Event analysis */

ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Check that the set of genes with MISO SE are expressed in our data */

*MISO genes;

data miso_se;
   set event.miso_se2gene;
   keep ens_gene_id;
run;

data ens2refseq;
   set event.ensembl2refseq_gene_id;
   keep ens_gene_id gene_id;
run;

data genes2keep;
   set event.flag_gene_expressed;
   where flag_gene_expressed=1;
   keep gene_id;
run;

proc sort data=ens2refseq nodup;
   by gene_id;
proc sort data=genes2keep nodup;
   by gene_id;
run;

data ens2refseq2;
  merge ens2refseq (in=in1) genes2keep (in=in2);
  by gene_id;
  if in1 and in2;
run; *14593 of 20875 with ENSG match;



proc sort data=miso_se nodup;
   by ens_gene_id;
proc sort data=ens2refseq2 nodup;
   by ens_gene_id gene_id;
run;

data miso_se2refseq;
   merge miso_se (in=in1) ens2refseq2 (in=in2);
   by ens_gene_id;
   if in1 and in2 then flag_has_refseq=1;
   else if in1 then flag_has_refseq=0;
   if in1 then output;
run;

proc freq data=miso_se2refseq;
  tables flag_has_refseq;
run;

/*
                                           Cumulative    Cumulative
lag_has_refseq    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0        1123       44.25          1123        44.25
             1        1415       55.75          2538       100.00

+ 13178 genes expressed without MISO

*/

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
         1 |   1123 |   1415 |   2538
           |  44.25 |  55.75 | 100.00
           |  44.25 |  55.75 |
           | 100.00 | 100.00 |
  ---------+--------+--------+
  Total        1123     1415     2538
              44.25    55.75   100.00

 flag_has_miso_se
           flag_has_refseq

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       1|  Total
 ---------+--------+
        1 |   1398 |   1398
          | 100.00 | 100.00
          | 100.00 |
          | 100.00 |
 ---------+--------+
 Total        1398     1398
            100.00   100.00

  flag_has_miso_se
            flag_has_multiple_exons

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         1 |     17 |   1398 |   1415
           |   1.20 |  98.80 | 100.00
           |   1.20 |  98.80 |
           | 100.00 | 100.00 |
  ---------+--------+--------+
  Total          17     1398     1415
               1.20    98.80   100.00


*/

data single_fusion;
   set miso2es;
   where flag_has_multiple_exons=0 and flag_has_refseq=1 and flag_has_miso_se=1;
run; *17 genes that have only a single exonic region, but have a MISO skipped exon;


/* make permenant */

data event.miso_refseq_exonskip_compare;
   set miso2es;
run;

 
