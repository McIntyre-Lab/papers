ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/*
GENE-LEVEL:
- Flag if gene has unannotated features detected
- Flag if gene has novel PB transcripts
- Check overlap
- Done!
*/

/* For each expressed gene : (global first)
 - flag if it has features with BLAST hits
 - flag if it has unannotated features
 - flag if it has unannotated features with BLAST hits
 - flag if it has PB transcripts
 - flag if it has novel PB transcripts
 - then look at overlaps! */

data feat2exp_genes;
  set event.feature2xs2gene_exp_only;
  drop transcript_id;
run;

data feat_data;
  set event.features_w_annotations;
  where flag_feature_on=1 and feature_type ne "fusion";
run;

proc sort data=feat2exp_genes nodup;
  by feature_id;
proc sort data=feat_data;
  by feature_id;
run;

data feat2exp_genes2;
   merge feat2exp_genes (in=in1) feat_data (in=in2);
   by feature_id;
   if in1 and in2;
run;

/* Genes with unannotated features detected */
data genes_w_unannot;
   set feat2exp_genes2;
   where flag_junction_annotated=0;
   keep gene_id;
run;

proc sort data=genes_w_unannot nodup;
  by gene_id;
run;

/* Genes with PacBio BLAST hits */

data frags_w_hits;
   set event.blast_pacbio_fragments_in_refseq;
   if flag_frag_redundant=1 then delete;
   if flag_annot_hit=1 or flag_unannot_hit=1 or flag_novel_hit=1;
   keep fragment_id;
   rename fragment_id=feature_id;
run;

data events_w_hits;
   set event.blast_pacbio_events_in_refseq;
   if flag_event_redundant=1 then delete;
   rename event_id=feature_id;
run;

data features_w_hits;
   set events_w_hits frags_w_hits;
run;

proc sort data=features_w_hits nodup;
   by feature_id;
proc sort data=feat2exp_genes2;
   by feature_id;
run;

data gene2feat_w_hits;
  merge features_w_hits (in=in1) feat2exp_genes2 (in=in2);
  by feature_id;
  if in1 and in2;
run;

data genes_w_hits;
  set gene2feat_W_hits;
  keep gene_id;
run;

proc sort data=genes_w_hits nodup;
  by gene_id;
run;

/* Genes with unannotated features with PacBio BLAST hits */
data unannot_feat_w_hits;
   set gene2feat_w_hits;
   if flag_junction_annotated=0;
   keep gene_id;
run;

proc sort data=unannot_feat_w_hits nodup;
  by gene_id;
run;

/* Genes with PB */

data gene2pb;
  set event.pacbio2refseq_gene;
  keep gene_id pacbio_gene_id;
run;

proc sort data=gene2pb nodup;
  by gene_id;
run;

/* Genes with novel PB */

data pb_gene2xs;
  set event.pacbio_transcripts;
  where pb_status="Novel";
  keep pacbio_gene_id;
run;

proc sort data=pb_gene2xs nodup;
  by pacbio_gene_id;
proc sort data=gene2pb ;
  by pacbio_gene_id;
run;

data pb_gene_w_novel;
  merge gene2pb (in=in1) pb_gene2xs (in=in2);
  by pacbio_gene_id;
  if in1 and in2;
  keep gene_id;
run;

/* Merge all together, then count */

proc sort data=genes_w_unannot nodup; *genes with unannot features;
  by gene_id;
proc sort data=genes_w_hits nodup; *genes with PacBio BLAST hits;
  by gene_id;
proc sort data=unannot_feat_w_hits nodup; *genes with unannot features with PacBio BLAST hits;
  by gene_id;
proc sort data=gene2pb nodup; *genes with PacBio transcripts;
  by gene_id;
proc sort data=pb_gene_w_novel nodup; *genes with novel PacBio transcripts;
  by gene_id;
run;


data all_gene_flags;
  merge genes_w_unannot (in=in1) genes_w_hits (in=in2) unannot_feat_w_hits (in=in3)
        gene2pb (in=in4) pb_gene_w_novel (in=in5);
  by gene_id;
  if in1 then flag_gene_has_unannot_feat=1; else flag_gene_has_unannot_feat=0;
  if in2 then flag_gene_has_pb_hits=1; else flag_gene_has_pb_hits=0; 
  if in3 then flag_gene_has_unannot_pb_hits=1; else flag_gene_has_unannot_pb_hits=0;
  if in4 then flag_gene_has_pb_xscripts=1; else flag_gene_has_pb_xscripts=0;
  if in5 then flag_gene_has_novel_pb=1; else flag_gene_has_novel_pb=0;
run;

/* Make permenant */

data event.genes_w_unannot_and_pb_hits;
   set all_gene_flags;
run;


proc freq data=event.genes_w_unannot_and_pb_hits noprint;
  tables flag_gene_has_unannot_feat*flag_gene_has_pb_hits*flag_gene_has_unannot_pb_hits*
         flag_gene_has_pb_xscripts*flag_gene_has_novel_pb / out=gene_blast_hits;
run;

proc print data=gene_blast_hits;
run;

/*
 flag_gene_     flag_gene_     flag_gene_     flag_gene_    flag_gene_
has_unannot_      has_pb_     has_unannot_      has_pb_     has_novel_
    feat           hits          pb_hits       xscripts         pb        COUNT

      0              0              0              1             0           51
      0              0              0              1             1           25
      1              0              0              0             0        19694
      1              0              0              1             0          448
      1              0              0              1             1           39
      1              1              0              0             0          194
      1              1              0              1             0         1724
      1              1              0              1             1         1084
      1              1              1              0             0          470
      1              1              1              1             0         1901
      1              1              1              1             1         1564

*/

/* I want to know:
Unannot vs unannot w/ hit
Unannot vs PB
Hit vs PB
Unannot vs Novel PB
Unannot w/ hit vs PB
Unannot w/ hit vs Novel PB
*/


proc freq data=event.genes_w_unannot_and_pb_hits;
  tables flag_gene_has_unannot_feat*flag_gene_has_unannot_pb_hits
         flag_gene_has_unannot_feat*flag_gene_has_pb_xscripts
         flag_gene_has_pb_hits*flag_gene_has_pb_xscripts
         flag_gene_has_unannot_feat*flag_gene_has_novel_pb
         flag_gene_has_unannot_pb_hits*flag_gene_has_pb_xscripts
         flag_gene_has_unannot_pb_hits*flag_gene_has_novel_pb;
run;



/*  Unannot vs unannot w/ hit:

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |     76 |      0 |     76
           |   0.28 |   0.00 |   0.28
           | 100.00 |   0.00 |
           |   0.33 |   0.00 |
  ---------+--------+--------+
         1 |  23183 |   3935 |  27118
           |  85.25 |  14.47 |  99.72
           |  85.49 |  14.51 |
           |  99.67 | 100.00 |
  ---------+--------+--------+
  Total       23259     3935    27194
              85.53    14.47   100.00

Unannot vs PB:

flag_gene_has_unannot_feat
          flag_gene_has_pb_xscripts

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |      0 |     76 |     76
         |   0.00 |   0.28 |   0.28
         |   0.00 | 100.00 |
         |   0.00 |   1.11 |
---------+--------+--------+
       1 |  20358 |   6760 |  27118
         |  74.86 |  24.86 |  99.72
         |  75.07 |  24.93 |
         | 100.00 |  98.89 |
---------+--------+--------+
Total       20358     6836    27194
            74.86    25.14   100.00

Hit vs PB : 
   flag_gene_has_pb_hits
             flag_gene_has_pb_xscripts

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  19694 |    563 |  20257
            |  72.42 |   2.07 |  74.49
            |  97.22 |   2.78 |
            |  96.74 |   8.24 |
   ---------+--------+--------+
          1 |    664 |   6273 |   6937
            |   2.44 |  23.07 |  25.51
            |   9.57 |  90.43 |
            |   3.26 |  91.76 |
   ---------+--------+--------+
   Total       20358     6836    27194
               74.86    25.14   100.00


Unannot vs Novel PB : 

              The SAS System           08:50 Th

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |     51 |     25 |     76
           |   0.19 |   0.09 |   0.28
           |  67.11 |  32.89 |
           |   0.21 |   0.92 |
  ---------+--------+--------+
         1 |  24431 |   2687 |  27118
           |  89.84 |   9.88 |  99.72
           |  90.09 |   9.91 |
           |  99.79 |  99.08 |
  ---------+--------+--------+
  Total       24482     2712    27194
              90.03     9.97   100.00

             The SAS System           08:50 T

Unannot w/ hit vs PB :
 flag_gene_has_unannot_pb_hits
           flag_gene_has_pb_xscripts

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  19888 |   3371 |  23259
          |  73.13 |  12.40 |  85.53
          |  85.51 |  14.49 |
          |  97.69 |  49.31 |
 ---------+--------+--------+
        1 |    470 |   3465 |   3935
          |   1.73 |  12.74 |  14.47
          |  11.94 |  88.06 |
          |   2.31 |  50.69 |
 ---------+--------+--------+
 Total       20358     6836    27194
             74.86    25.14   100.00


Unannot w/ hit vs Novel PB :
flag_gene_has_unannot_pb_hits
            flag_gene_has_novel_pb

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  22111 |   1148 |  23259
           |  81.31 |   4.22 |  85.53
           |  95.06 |   4.94 |
           |  90.32 |  42.33 |
  ---------+--------+--------+
         1 |   2371 |   1564 |   3935
           |   8.72 |   5.75 |  14.47
           |  60.25 |  39.75 |
           |   9.68 |  57.67 |
  ---------+--------+--------+
  Total       24482     2712    27194
              90.03     9.97   100.00

About half. Possible that novel PB transcripts encompass more than we can capture using 
  generated events

*/  







