ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* What is the overlap between genes with DD exons and genes with DD PB isoforms? */

data pb2refseq_as;
   set event.pb2refseq_flag_gene_alt_spliced;
   keep gene_id flag_pb2rs_single_iso_gene flag_pb2rs_gene_as flag_pb2rs_gene_specific;
run;

data event_as;
  set event.genes_all_events_on_by_cell;
  keep gene_id flag_single_fusion_gene flag_gene_fusions_as flag_gene_fusions_specific flag_gene_fusions_off;
run;

proc sort data=pb2refseq_as;
  by gene_id;
proc sort data=event_as;
  by gene_id;
run;

data pb2event_as;
  merge pb2refseq_as (in=in1) event_as (in=in2);
  by gene_id;
  if in1 and in2;
run; *4719 genes in common;

proc freq data=pb2event_as;
  where flag_pb2rs_single_iso_gene=0 and flag_pb2rs_gene_specific=0
    and flag_single_fusion_gene=0 and flag_gene_fusions_specific=0 and flag_gene_fusions_off=0;
  tables flag_gene_fusions_as*flag_pb2rs_gene_as;
run;

/*
  flag_gene_fusions_as
            flag_pb2rs_gene_as

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |   1365 |    360 |   1725
           |  62.99 |  16.61 |  79.60
           |  79.13 |  20.87 |
           |  80.63 |  75.95 |
  ---------+--------+--------+
         1 |    328 |    114 |    442
           |  15.14 |   5.26 |  20.40
           |  74.21 |  25.79 |
           |  19.37 |  24.05 |
  ---------+--------+--------+
  Total        1693      474     2167
              78.13    21.87   100.00
*/


proc freq data=pb2event_as noprint;
  tables flag_pb2rs_single_iso_gene*flag_pb2rs_gene_specific*flag_single_fusion_gene*
         flag_gene_fusions_specific*flag_gene_fusions_off*flag_gene_fusions_as*flag_pb2rs_gene_as / out=as_count;
proc print data=as_count;
run;

/*
 flag_pb2rs_                            flag_gene_ flag_gene_             flag_
 single_iso_  flag_pb2rs_  flag_single_  fusions_   fusions_  flag_gene_  pb2rs_
     gene    gene_specific  fusion_gene  specific      off    fusions_as gene_as COUNT

      0            0             0           0          0          0        0     1365
      0            0             0           0          0          0        1      360
      0            0             0           0          0          1        0      328
      0            0             0           0          0          1        1      114
      0            0             1           0          0          0        0        4
      0            1             0           0          0          0        0       14
      0            1             0           0          0          1        0       27
      1            0             0           0          0          0        0     1742
      1            0             0           0          0          1        0      579
      1            0             1           0          0          0        0       70
      1            1             0           0          0          0        0       23
      1            1             0           0          0          1        0       89
      1            1             1           0          0          0        0        4
*/



data pb2refseq_as;
   set event.pb2refseq_flag_gene_alt_spliced;
   keep gene_id flag_pb2rs_single_iso_gene flag_pb2rs_gene_as flag_pb2rs_gene_specific;
run;

data event_as;
  set event.genes_all_events_on_by_cell;
  keep gene_id flag_single_fragment_gene flag_gene_fragments_as flag_gene_fragments_specific flag_gene_fragments_off;
run;

proc sort data=pb2refseq_as;
  by gene_id;
proc sort data=event_as;
  by gene_id;
run;

data pb2event_as;
  merge pb2refseq_as (in=in1) event_as (in=in2);
  by gene_id;
  if in1 and in2;
run; *4719 genes in common;

proc freq data=pb2event_as;
  where flag_pb2rs_single_iso_gene=0 and flag_pb2rs_gene_specific=0
    and flag_single_fragment_gene=0 and flag_gene_fragments_specific=0 and flag_gene_fragments_off=0;
  tables flag_gene_fragments_as*flag_pb2rs_gene_as;
run;

/*
flag_gene_fragments_as
          flag_pb2rs_gene_as

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   1223 |    305 |   1528
         |  56.44 |  14.07 |  70.51
         |  80.04 |  19.96 |
         |  72.24 |  64.35 |
---------+--------+--------+
       1 |    470 |    169 |    639
         |  21.69 |   7.80 |  29.49
         |  73.55 |  26.45 |
         |  27.76 |  35.65 |
---------+--------+--------+
Total        1693      474     2167
            78.13    21.87   100.00
*/




/* What is the overlap between genes with DD exons and genes with DD PB isoforms?
   Let's also remove RefSeq genes that have more than 1 PacBio gene assigned to them (as these could be ambiguous)
*/

data pb2refseq_as;
   set event.pb2refseq_flag_gene_alt_spliced;
   if num_pb_genes > 1 then delete;
   keep gene_id flag_pb_single_iso_gene flag_pb_gene_as flag_pb_gene_specific;
run;

data event_as;
  set event.genes_all_events_on_by_cell;
  keep gene_id flag_single_fusion_gene flag_gene_fusions_as flag_gene_fusions_specific flag_gene_fusions_off;
run;

proc sort data=pb2refseq_as;
  by gene_id;
proc sort data=event_as;
  by gene_id;
run;

data pb2event_as;
  merge pb2refseq_as (in=in1) event_as (in=in2);
  by gene_id;
  if in1 and in2;
run; *4561 genes in common;

proc freq data=pb2event_as;
  where flag_pb_single_iso_gene=0 and flag_pb_gene_specific=0
    and flag_single_fusion_gene=0 and flag_gene_fusions_specific=0 and flag_gene_fusions_off=0;
  tables flag_gene_fusions_as*flag_pb_gene_as;
run;

/*
flag_gene_fusions_as
          flag_pb_gene_as

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   1288 |    326 |   1614
         |  64.02 |  16.20 |  80.22
         |  79.80 |  20.20 |
         |  81.36 |  75.99 |
---------+--------+--------+
       1 |    295 |    103 |    398
         |  14.66 |   5.12 |  19.78
         |  74.12 |  25.88 |
         |  18.64 |  24.01 |
---------+--------+--------+
Total        1583      429     2012
            78.68    21.32   100.00


Note: this isn't taking into consideration isoforms that might utilize the same set of exons, so I don't
think the comparison is particularly good
*/


proc freq data=pb2event_as noprint;
  tables flag_pb_single_iso_gene*flag_pb_gene_specific*flag_single_fusion_gene*
         flag_gene_fusions_specific*flag_gene_fusions_off*flag_gene_fusions_as*flag_pb_gene_as / out=as_count;
proc print data=as_count;
run;

/*
flag_pb_                            flag_gene_ flag_gene_
 single_ flag_pb_gene_ flag_single_  fusions_   fusions_  flag_gene_ flag_pb_
iso_gene    specific    fusion_gene  specific      off    fusions_as  gene_as COUNT PERCENT

    0          0             0           0          0          0         0     1288 28.2394
    0          0             0           0          0          0         1      326  7.1476
    0          0             0           0          0          1         0      295  6.4679
    0          0             0           0          0          1         1      103  2.2583
    0          0             1           0          0          0         0        1  0.0219
    0          1             0           0          0          0         0       14  0.3070
    0          1             0           0          0          1         0       27  0.5920
    1          0             0           0          0          0         0     1742 38.1934
    1          0             0           0          0          1         0      579 12.6946
    1          0             1           0          0          0         0       70  1.5348
    1          1             0           0          0          0         0       23  0.5043
    1          1             0           0          0          1         0       89  1.9513
    1          1             1           0          0          0         0        4  0.0877
*/



data pb2refseq_as;
   set event.pb2refseq_flag_gene_alt_spliced;
   if num_pb_genes > 1 then delete;
   keep gene_id flag_pb_single_iso_gene flag_pb_gene_as flag_pb_gene_specific;
run;


data event_as;
  set event.genes_all_events_on_by_cell;
  keep gene_id flag_single_fragment_gene flag_gene_fragments_as flag_gene_fragments_specific flag_gene_fragments_off;
run;

proc sort data=pb2refseq_as;
  by gene_id;
proc sort data=event_as;
  by gene_id;
run;

data pb2event_as;
  merge pb2refseq_as (in=in1) event_as (in=in2);
  by gene_id;
  if in1 and in2;
run; *4561 genes in common;

proc freq data=pb2event_as;
  where flag_pb_single_iso_gene=0 and flag_pb_gene_specific=0
    and flag_single_fragment_gene=0 and flag_gene_fragments_specific=0 and flag_gene_fragments_off=0;
  tables flag_gene_fragments_as*flag_pb_gene_as;
run;

/*
flag_gene_fragments_as
          flag_pb_gene_as

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   1161 |    279 |   1440
         |  57.70 |  13.87 |  71.57
         |  80.63 |  19.38 |
         |  73.34 |  65.03 |
---------+--------+--------+
       1 |    422 |    150 |    572
         |  20.97 |   7.46 |  28.43
         |  73.78 |  26.22 |
         |  26.66 |  34.97 |
---------+--------+--------+
Total        1583      429     2012
            78.68    21.32   100.00

Use for an example
A= E+ and P+ = 150
B= E- and P+ = 279
C= E+ and P- = 422
D= E- and P- = 1161

FP = C/(C+D) = 422/(422+1161) = 17%
FN = B/(A+B) = 279/(279+150) = 65%

*/


