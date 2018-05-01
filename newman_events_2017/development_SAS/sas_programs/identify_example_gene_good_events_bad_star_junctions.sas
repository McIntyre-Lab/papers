/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Identify genes that Event Analysis performs well for, but poorly for STAR
   (1) Only consider expressed genes without multigene exonic regions
   (2) Count the number of detected annotated junctions per gene
   (3) Count the number of STAR detected annotated junctions per gene
*/

data gene2keep;
  set event.feature2xs2gene_exp_only_nomulti;
  keep gene_id;
run;

data annot;
   set evspl.splicing_events_annot_refseq;
   where flag_junction_annotated=1 or num_transcripts > 0;
   keep gene_id event_id;
run;

proc sort data=gene2keep nodup;
   by gene_id;
proc sort data=annot;
   by gene_id event_id;
run;

data annot2keep;
   merge gene2keep (in=in1) annot (in=in2);
   by gene_id;
   if in1 and in2;
run;


data junc2event_annot;
   set eventloc.unique_junction2event_mm10;
   keep event_id junction_id;
run;

proc sort data=annot2keep;
  by event_id;
proc sort data=junc2event_annot nodup;
  by event_id junction_id;
run;

data annot2keep2;
  merge annot2keep (in=in1) junc2event_annot (in=in2);
  by event_id;
  if in1 and in2;
run;

/* Junction table */

data junc_table;
  set event.event2star2pacbio_junc_table;
  where flag_in_catalog=1;
  keep junction_id flag_events_NSC_apn_gt0 flag_events_NSC_apn_ge5
                   flag_star_NSC_apn_gt0 flag_star_NSC_apn_ge5;
run;

proc sort data=annot2keep2;
  by junction_id;
proc sort data=junc_table;
  by junction_id;
run;

data annot2keep_support;
  merge annot2keep2 (in=in1) junc_table (in=in2);
  by junction_id;
  if in1 and in2;
run;

/* Per gene, calc junctions with support */

proc sort data=annot2keep_support;
  by gene_id junction_id;
proc means data=annot2keep_support noprint;
  by gene_id;
  var flag_events_NSC_apn_gt0 flag_events_NSC_apn_ge5
                   flag_star_NSC_apn_gt0 flag_star_NSC_apn_ge5;
  output out=num_junc_dtct_by_gene
         sum(flag_events_NSC_apn_gt0)=num_events_junc_support_apn0
         sum(flag_events_NSC_apn_ge5)=num_events_junc_support_apn5
         sum(flag_star_NSC_apn_gt0)=num_star_junc_support_apn0
         sum(flag_star_NSC_apn_ge5)=num_star_junc_support_apn5;
run;

data perc_junc_on_by_gene;
  set num_junc_dtct_by_gene;
  perc_events_junc_support_apn0=num_events_junc_support_apn0/_FREQ_;
  perc_events_junc_support_apn5=num_events_junc_support_apn5/_FREQ_;
  perc_star_junc_support_apn0=num_star_junc_support_apn0/_FREQ_;
  perc_star_junc_support_apn5=num_star_junc_support_apn5/_FREQ_;
  drop _TYPE_;
  rename _FREQ_=total_junctions;
run;

/* Flag genes with all catalog junctions supported,
   greater than 50% of STAR junctions supported */

data flag_genes;
  set perc_junc_on_by_gene;
  if total_junctions > 5 then flag_num_junctions_gt5=1;
  else flag_num_junctions_gt5=0;
  if perc_events_junc_support_apn0 = 1 then flag_events_all_junc_apn0=1;
  else flag_events_all_junc_apn0=0;
  if perc_star_junc_support_apn0 ge 0.5 then flag_star_50perc_junc_apn0=1;
  else flag_star_50perc_junc_apn0=0;

  if perc_events_junc_support_apn5 = 1 then flag_events_all_junc_apn5=1;
  else flag_events_all_junc_apn5=0;
  if perc_star_junc_support_apn5 ge 0.5 then flag_star_50perc_junc_apn5=1;
  else flag_star_50perc_junc_apn5=0;
run;

proc freq data=flag_genes noprint;
  tables flag_num_junctions_gt5*flag_events_all_junc_apn0*flag_star_50perc_junc_apn0
        / out=junc_summary_apn0;

  tables flag_num_junctions_gt5*flag_events_all_junc_apn5*flag_star_50perc_junc_apn5
        / out=junc_summary_apn5;

  tables flag_num_junctions_gt5*flag_events_all_junc_apn5*flag_star_50perc_junc_apn0
        / out=junc_summary_apn5_star_apn0;
run;

proc print data=junc_summary_apn0;
run;
/*
   flag_num_    flag_events_    flag_star_
  junctions_      all_junc_       50perc_
      gt5           apn0         junc_apn0    COUNT

       0              0              0         7476
       0              0              1          337
       0              1              0          304
       0              1              1          933
       1              0              0         6565
       1              0              1         2494
       1              1              0           55
       1              1              1          927

55 possible genes here
*/

proc print data=junc_summary_apn5;
run;

/*
  flag_num_    flag_events_    flag_star_
 junctions_      all_junc_       50perc_
     gt5           apn5         junc_apn5    COUNT

      0              0              0         8368
      0              0              1            7
      0              1              0          638
      0              1              1           37
      1              0              0         9636
      1              0              1            6
      1              1              0          386
      1              1              1           13

386 possible genes here
*/

proc print data=junc_summary_apn5_star_apn0;
run;
/*
 flag_num_    flag_events_    flag_star_
junctions_      all_junc_       50perc_
    gt5           apn5         junc_apn0    COUNT

     0              0              0         7665
     0              0              1          710
     0              1              0          115
     0              1              1          560
     1              0              0         6590
     1              0              1         3052
     1              1              0           30
     1              1              1          369

30 possible genes here
*/

/* Pull out the 30 genes that have good junction support at APN5 in events,
   but poor junction support at APN0 using STAR */

data genes4fig;
  set flag_genes;
  where flaG_num_junctions_gt5=1 and flag_events_all_junc_apn5=1 and flag_star_50perc_junc_apn0=0;
run;

/* Repeat above, but for NIC junctions. I want genes with at least 3 NIC junctions detected */



data gene2keep;
  set event.feature2xs2gene_exp_only_nomulti;
  keep gene_id;
run;

data nic;
   set evspl.splicing_events_annot_refseq;
   where flag_junction_annotated=0 and num_transcripts = 0 and flag_intron_retention=0;
   keep gene_id event_id;
run;

proc sort data=gene2keep nodup;
   by gene_id;
proc sort data=nic;
   by gene_id event_id;
run;

data nic2keep;
   merge gene2keep (in=in1) nic (in=in2);
   by gene_id;
   if in1 and in2;
run;

data junc2event_nic;
   set eventloc.unique_junction2event_mm10;
   keep event_id junction_id;
run;

proc sort data=nic2keep;
  by event_id;
proc sort data=junc2event_nic nodup;
  by event_id junction_id;
run;

data nic2keep2;
  merge nic2keep (in=in1) junc2event_nic (in=in2);
  by event_id;
  if in1 and in2;
run;

/* Junction table */

data junc_table;
  set event.event2star2pacbio_junc_table;
  where flag_in_catalog=1 and junction_type ? "NIC";
  keep junction_id flag_events_NSC_apn_gt0 flag_events_NSC_apn_ge5
                   flag_star_NSC_apn_gt0 flag_star_NSC_apn_ge5;
run;

proc sort data=nic2keep2;
  by junction_id;
proc sort data=junc_table;
  by junction_id;
run;

data nic2keep_support;
  merge nic2keep2 (in=in1) junc_table (in=in2);
  by junction_id;
  if in1 and in2;
run;

/* Per gene, calc junctions with support */

proc sort data=nic2keep_support;
  by gene_id junction_id;
proc means data=nic2keep_support noprint;
  by gene_id;
  var flag_events_NSC_apn_gt0 flag_events_NSC_apn_ge5
                   flag_star_NSC_apn_gt0 flag_star_NSC_apn_ge5;
  output out=num_nic_dtct_by_gene
         sum(flag_events_NSC_apn_gt0)=num_events_nic_support_apn0
         sum(flag_events_NSC_apn_ge5)=num_events_nic_support_apn5
         sum(flag_star_NSC_apn_gt0)=num_star_nic_support_apn0
         sum(flag_star_NSC_apn_ge5)=num_star_nic_support_apn5;
run;

data flag_nic_on_by_gene;
  set num_nic_dtct_by_gene;
  if num_events_nic_support_apn0 ge 3 then flag_events_nic_gt3_apn0=1;
  else flag_events_nic_gt3_apn0=0;

  if num_events_nic_support_apn5 ge 3 then flag_events_nic_gt3_apn5=1;
  else flag_events_nic_gt3_apn5=0;

  if num_star_nic_support_apn0 ge 1 then flag_star_nic_gt1_apn0=1;
  else flag_star_nic_gt1_apn0=0;

  if num_star_nic_support_apn5 ge 1 then flag_star_nic_gt1_apn5=1;
  else flag_star_nic_gt1_apn5=0;

  drop _TYPE_ _FREQ_;
run;

proc freq data=flag_nic_on_by_gene;
  tables flag_events_nic_gt3_apn5*flag_star_nic_gt1_apn0;
run;

/*
  flag_events_nic_gt3_apn5
            flag_star_nic_gt1_apn0

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  16035 |    454 |  16489
           |  97.24 |   2.75 |  99.99
           |  97.25 |   2.75 |
           | 100.00 |  99.78 |
  ---------+--------+--------+
         1 |      0 |      1 |      1
           |   0.00 |   0.01 |   0.01
           |   0.00 | 100.00 |
           |   0.00 |   0.22 |
  ---------+--------+--------+
  Total       16035      455    16490
              97.24     2.76   100.00
*/


proc freq data=flag_nic_on_by_gene;
  tables flag_events_nic_gt3_apn0*flag_star_nic_gt1_apn0;
run;

/*
   flag_events_nic_gt3_apn0
             flag_star_nic_gt1_apn0

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |  15981 |    395 |  16376
            |  96.91 |   2.40 |  99.31
            |  97.59 |   2.41 |
            |  99.66 |  86.81 |
   ---------+--------+--------+
          1 |     54 |     60 |    114
            |   0.33 |   0.36 |   0.69
            |  47.37 |  52.63 |
            |   0.34 |  13.19 |
   ---------+--------+--------+
   Total       16035      455    16490
               97.24     2.76   100.00

Pick one of these 54 genes.
*/

data nic_genes;
  set flag_nic_on_by_gene;
  if flag_events_nic_gt3_apn0=1 and flag_star_nic_gt1_apn0=0;
  keep gene_id;
run;

proc sort data=nic_genes;
  by gene_id;
proc sort data=genes4fig;
  by gene_id;
run;

data genes4gfig_2;
  merge genes4fig (in=in1) nic_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* 1 gene remaining! -- geneID=27367, 11 junctions total, only 1 detected with STAR, only 2 transcripts */



data nic_genes;
  set flag_nic_on_by_gene;
  if flag_events_nic_gt3_apn0=1 ;
  keep gene_id;
run;

proc sort data=nic_genes;
  by gene_id;
proc sort data=genes4fig;
  by gene_id;
run;

data genes4gfig_3;
  merge genes4fig (in=in1) nic_genes (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* 3 genes:
19167	12 junctions, 4 STAR, 3-4 transcripts
27367	11 junctions, 1 STAR, 2 transcripts
53600	6 junctions, 0 STAR, 1 transcript

Pick the first two, let LMM decide

*/
