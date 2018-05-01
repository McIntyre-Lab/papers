ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Check that the set of genes with MISO SE are the same set with SE in Event analysis 
   I want to flag genes that have events in both, so we can focus on just this list*/



*MISO genes;
data miso_se;
   set event.miso_se2gene;
   keep ens_gene_id;
run;

data ens2refseq;
   set event.ensembl2refseq_gene_id;
   keep ens_gene_id gene_id;
run;

proc sort data=miso_se nodup;
   by ens_gene_id;
proc sort data=ens2refseq nodup;
   by ens_gene_id gene_id;
run;

data miso_se2refseq;
   merge miso_se (in=in1) ens2refseq (in=in2);
   by ens_gene_id;
   if in1 and in2 then flag_has_refseq=1;
   else if in1 then flag_has_refseq=0;
   if in1 then output;
run;

data miso_events;
   set event.miso_se2gene;
run;

data miso_genes;
   set miso_se2refseq;
   where flag_has_refseq=1;
   keep ens_gene_id;
run;


proc sort data=miso_genes nodup;
   by ens_gene_id;
proc sort data=miso_events;
   by ens_gene_id;
run;

data miso_events_gene_in_refseq;
  merge miso_events (in=in1) miso_genes (in=in2);
  by ens_gene_id;
  if in1 and in2;
run;
*3187 genes with Refseq match;


*Exonskipping genes;

data exonskip;
  set evspl.splicing_events_annot_refseq;
  where flag_exonskip=1;
  keep gene_id;
run;

proc sort data=miso_se2refseq nodup;
   by gene_id;
proc sort data=exonskip nodup;
   by gene_id;
run;

data miso2es;
  merge miso_se2refseq (in=in1) exonskip (in=in2);
   by gene_id;
  if in1 then flag_has_miso_se=1; else flag_has_miso_se=0;
  if in2 then flag_has_exonskip_event=1; else flag_has_exonskip_event=0;
run;

proc freq data=miso2es;
  tables flag_has_miso_se*flag_has_exonskip_event;
run;

/*
 flag_has_miso_se
           flag_has_exonskip_event

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |  26800 |  26800
          |   0.00 |  91.30 |  91.30
          |   0.00 | 100.00 |
          |   0.00 |  93.01 |
 ---------+--------+--------+
        1 |    541 |   2014 |   2555
          |   1.84 |   6.86 |   8.70
          |  21.17 |  78.83 |
          | 100.00 |   6.99 |
 ---------+--------+--------+
 Total         541    28814    29355
              1.84    98.16   100.00
*/

proc freq data=miso2es;
  where flag_has_refseq=1;
  tables flag_has_miso_se*flag_has_exonskip_event;
run;

/*

 flag_has_miso_se
           flag_has_exonskip_event

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        1 |     87 |   2014 |   2101
          |   4.14 |  95.86 | 100.00
          |   4.14 |  95.86 |
          | 100.00 | 100.00 |
 ---------+--------+--------+
 Total          87     2014     2101
              4.14    95.86   100.00

*/

data no_exonskip;
   set miso2es;
   where flag_has_exonskip_event=0 and flag_has_refseq=1 and flag_has_miso_se=1;
run; *87 genes that don't have an exonskipping junctions, but have a MISO skipped exon;


/* make permenant */

data event.miso_refseq_exonskip_compare;
   set miso2es;
run;

 data event.miso_events_gene_in_refseq;
   set miso_events_gene_in_refseq;
run;

