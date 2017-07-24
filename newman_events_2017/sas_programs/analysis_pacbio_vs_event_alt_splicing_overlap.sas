ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';

/* Compare PB genes with differentially detected isoforms against 
   Event Analysis genes that are alternatively spliced */

* PB2Refseq genes, flagged AS;

data pacbio_as;
   set event.pacbio2refseq_flag_diff_xs;
run;

* Event Analysis genes, flagged AS;
data event_as;
   set event.flag_gene_alt_spliced;
run;

/* Merge */

proc sort data=pacbio_as;
   by gene_id;
proc sort data=event_as;
   by gene_id;
run;

data pacbio2event_as;
  merge pacbio_as (in=in1) event_as (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* For Event genes with Pacbio sequences, what is the % of genes alternatively spliced ? */

proc freq data=pacbio2event_as;
   where flag_gene_has_pb=1 and flag_pb_gene_multiple_xs=1;
   tables flag_pb_gene_as flag_gene_fusion_dd flag_gene_fragment_dd flag_gene_splicing_dd flag_gene_as;
run;

/*
                                              Cumulative    Cumulative
 flag_pb_gene_as    Frequency     Percent     Frequency      Percent
 --------------------------------------------------------------------
               0        1685       76.52          1685        76.52
               1         517       23.48          2202       100.00


      flag_gene_                             Cumulative    Cumulative
       fusion_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        1742       79.11          1742        79.11
               1         460       20.89          2202       100.00


      flag_gene_                             Cumulative    Cumulative
     fragment_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        1378       62.58          1378        62.58
               1         824       37.42          2202       100.00


      flag_gene_                             Cumulative    Cumulative
     splicing_dd    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0         106        4.81           106         4.81
               1        2096       95.19          2202       100.00


                                           Cumulative    Cumulative
  flag_gene_as    Frequency     Percent     Frequency      Percent
  -----------------------------------------------------------------
             0          95        4.31            95         4.31
             1        2107       95.69          2202       100.00
*/

proc freq data=pacbio2event_as;
   where flag_gene_has_pb=1 and flag_pb_gene_multiple_xs=1;
   tables flag_pb_gene_as*flag_gene_fusion_dd
          flag_pb_gene_as*flag_gene_fragment_dd
          flag_pb_gene_as*flag_gene_splicing_dd
          flag_pb_gene_as*flag_gene_as;
run;

/*
   flag_pb_gene_as
             flag_gene_fusion_dd

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |   1366 |    319 |   1685
            |  62.03 |  14.49 |  76.52
            |  81.07 |  18.93 |
            |  78.42 |  69.35 |
   ---------+--------+--------+
          1 |    376 |    141 |    517
            |  17.08 |   6.40 |  23.48
            |  72.73 |  27.27 |
            |  21.58 |  30.65 |
   ---------+--------+--------+
   Total        1742      460     2202
               79.11    20.89   100.00

 Table of flag_pb_gene_as by flag_gene_fragment_dd

        flag_pb_gene_as
                  flag_gene_fragment_dd

        Frequency|
        Percent  |
        Row Pct  |
        Col Pct  |       0|       1|  Total
        ---------+--------+--------+
               0 |   1100 |    585 |   1685
                 |  49.95 |  26.57 |  76.52
                 |  65.28 |  34.72 |
                 |  79.83 |  71.00 |
        ---------+--------+--------+
               1 |    278 |    239 |    517
                 |  12.62 |  10.85 |  23.48
                 |  53.77 |  46.23 |
                 |  20.17 |  29.00 |
        ---------+--------+--------+
        Total        1378      824     2202
                    62.58    37.42   100.00

  Table of flag_pb_gene_as by flag_gene_splicing_dd

         flag_pb_gene_as
                   flag_gene_splicing_dd

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 |     89 |   1596 |   1685
                  |   4.04 |  72.48 |  76.52
                  |   5.28 |  94.72 |
                  |  83.96 |  76.15 |
         ---------+--------+--------+
                1 |     17 |    500 |    517
                  |   0.77 |  22.71 |  23.48
                  |   3.29 |  96.71 |
                  |  16.04 |  23.85 |
         ---------+--------+--------+
         Total         106     2096     2202
                      4.81    95.19   100.00

  Table of flag_pb_gene_as by flag_gene_as

    flag_pb_gene_as     flag_gene_as

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 |     80 |   1605 |   1685
             |   3.63 |  72.89 |  76.52
             |   4.75 |  95.25 |
             |  84.21 |  76.17 |
    ---------+--------+--------+
           1 |     15 |    502 |    517
             |   0.68 |  22.80 |  23.48
             |   2.90 |  97.10 |
             |  15.79 |  23.83 |
    ---------+--------+--------+
    Total          95     2107     2202
                 4.31    95.69   100.00

*/

