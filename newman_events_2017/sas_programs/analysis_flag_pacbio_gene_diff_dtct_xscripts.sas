ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';

/* Identify PacBio genes that are expressing different isoforms in NPCs and OLDs
   This will be used for my "true positive" rate */

* PB transcript detection,add geneID;

data pb_dtct;
  length pacbio_gene_id $15.;
  set event.pacbio_xscripts_flag_cell_exp;
  pacbio_gene_id=catt("PB.",scan(pacbio_id,2,"."));
run;

data flag_pb_diff;
   set pb_dtct;
   if sum(flag_npc_expressed,flag_old_expressed) = 1 then flag_pb_xscript_dd=1;
   else flag_pb_xscript_dd=0;
run;

proc freq data=flag_pb_diff;
  tables flag_pb_xscript_dd;
run;

/*

flag_pb_xscript_                             Cumulative    Cumulative
              dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       13862       91.53         15145       100.00
               1        1283        8.47          1283         8.47
*/

proc sort data=flag_pb_diff;
   by pacbio_gene_id;
proc means data=flag_pb_diff noprint;
   by pacbio_gene_id;
   var flag_pb_xscript_dd;
   output out=flag_pb_gene_diff max(flag_pb_xscript_dd)=flag_pb_gene_as;
run;

data flag_pb_gene_diff2;
   set flag_pb_gene_diff;
   if _FREQ_ = 1 then flag_pb_gene_multiple_xs=0;
   else flag_pb_gene_multiple_xs=1;
run;


proc freq data=flag_pb_gene_diff2;
  tables flag_pb_gene_as flag_pb_gene_as*flag_pb_gene_multiple_xs;
run;

/*
                                              Cumulative    Cumulative
  flag_pb_gene_as    Frequency     Percent     Frequency      Percent
  --------------------------------------------------------------------
                0        6678       86.96          6678        86.96
                1        1001       13.04          7679       100.00

 flag_pb_gene_as
           flag_pb_gene_multiple_xs

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   4152 |   2526 |   6678
          |  54.07 |  32.89 |  86.96
          |  62.17 |  37.83 |
          |  94.36 |  77.04 |
 ---------+--------+--------+
        1 |    248 |    753 |   1001
          |   3.23 |   9.81 |  13.04
          |  24.78 |  75.22 |
          |   5.64 |  22.96 |
 ---------+--------+--------+
 Total        4400     3279     7679
             57.30    42.70   100.00

13% of genes, but when considering only genes that have more than 1 annotated isoform (!!)
this is 23% of genes
*/


/* PB2Refseq gene ID */

data pb2rs_gene;
   set event.pacbio_gene_to_refseq;
   keep gene_id pacbio_gene_id;
run;

proc sort data=pb2rs_gene nodup;
  by pacbio_gene_id gene_id;
proc sort data=flag_pb_gene_diff2;
   by pacbio_gene_id;
run;

data pb_diff_w_rs;
  merge flag_pb_gene_diff2 (in=in1) pb2rs_gene (in=in2);
   by pacbio_gene_id;
  if in1 and in2;
run;

proc freq data=pb_diff_w_rs;
  tables flag_pb_gene_as*flag_pb_gene_multiple_xs;
run;

/*
 flag_pb_gene_as
           flag_pb_gene_multiple_xs

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   4152 |   2558 |   6710
          |  53.70 |  33.08 |  86.78
          |  61.88 |  38.12 |
          |  94.36 |  76.77 |
 ---------+--------+--------+
        1 |    248 |    774 |   1022
          |   3.21 |  10.01 |  13.22
          |  24.27 |  75.73 |
          |   5.64 |  23.23 |
 ---------+--------+--------+
 Total        4400     3332     7732
             56.91    43.09   100.00

*/

/* Some PB genes map to multiple RefSeq genes. I am going to collapse PacBio gene AS to RefSeq
   I am also going to flag these genes in case we want to remove them later */

proc sort data=pb_diff_w_rs;
  by gene_id;
proc means data=pb_diff_w_rs noprint;
  by gene_id;
  var flag_pb_gene_as flag_pb_gene_multiple_xs;
  output out=pb_diff_rs_collapsed max=;
run;


/* Remove multigene */

data gene2keep;
   set event.flag_gene_expressed;
   where flag_gene_expressed=1;
   keep gene_id;
run;

proc sort data=gene2keep;
  by gene_id;
proc sort data=pb_diff_rs_collapsed;
  by gene_id;
run;

data pb_diff_w_rs_kept;
  merge gene2keep (in=in1) pb_diff_rs_collapsed (in=in2);
  by gene_id;
  if in1 and in2 then do;
    flag_gene_has_pb=1;
    if _FREQ_ gt 1 then flag_pb_multigene=1; else flag_pb_multigene=0;
    output; end;
  else if in1 then do;
    flag_gene_has_pb=0;
    output; end;
run;

proc freq data=pb_diff_w_rs_kept;
  where flag_gene_has_pb=1;
  tables flag_pb_gene_as*flag_pb_gene_multiple_xs;
run;

/*
  flag_pb_gene_as
            flag_pb_gene_multiple_xs

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   2436 |   1685 |   4121
           |  51.19 |  35.41 |  86.59
           |  59.11 |  40.89 |
           |  95.27 |  76.52 |
  ---------+--------+--------+
         1 |    121 |    517 |    638
           |   2.54 |  10.86 |  13.41
           |  18.97 |  81.03 |
           |   4.73 |  23.48 |
  ---------+--------+--------+
  Total        2557     2202     4759
              53.73    46.27   100.00

23% of genes with multiple detected isoforms

*/

proc freq data=pb_diff_w_rs_kept;
  where flag_gene_has_pb=1 and flag_pb_multigene=0;
  tables flag_pb_gene_as*flag_pb_gene_multiple_xs;
run;

/*
 flag_pb_gene_as
           flag_pb_gene_multiple_xs

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   2387 |   1618 |   4005
          |  51.96 |  35.22 |  87.18
          |  59.60 |  40.40 |
          |  95.37 |  77.38 |
 ---------+--------+--------+
        1 |    116 |    473 |    589
          |   2.53 |  10.30 |  12.82
          |  19.69 |  80.31 |
          |   4.63 |  22.62 |
 ---------+--------+--------+
 Total        2503     2091     4594
             54.48    45.52   100.00

*/

/* Make permenant */

data event.pacbio2refseq_flag_diff_xs;
  set pb_diff_w_rs_kept;
run;

