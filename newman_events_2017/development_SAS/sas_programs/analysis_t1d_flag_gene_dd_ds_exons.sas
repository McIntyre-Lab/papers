ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For the mouse data, with the set of genes that are expressed in both cell types:
(1) Count genes that have differentially detected exons (ie, it has exons that are only seen in NPCs
    and exons only seen in OPCs
(2) How many of these overlap with the number of genes tested for quant splice?

This will form the set of genes that are differentially or alternatively spliced (DAS) */

data flag_gene_dd;
   set event.t1d_flag_gene_on_by_cell;
   /* If there are more exons detected in either cell type than in both cell types together,
      then the gene has differentially detected exons */

    total_fus_per_gene=flag_fusion_cd4_on+flag_fusion_cd8_on+flag_fusion_cd19_on+
                       flag_fusion_cd4cd8_on+flag_fusion_cd4cd19_on+flag_fusion_cd8cd19_on+
                       flag_fusion_all_on;

    if total_fus_per_gene= 1 then flag_gene_monoexon=1;
    else flag_gene_monoexon=0;

   if flag_gene_cd4_on=1 and flag_gene_cd8_on=1 and flag_gene_cd19_on=1 then do;
      flag_gene_cell_specific=0;
      if flag_fusion_cd4_on > 0
         or flag_fusion_cd8_on > 0
         or flag_fusion_cd19_on > 0
         or flag_fusion_cd4cd8_on > 0
         or flag_fusion_cd4cd19_on > 0
         or flag_fusion_cd8cd19_on > 0
         then flag_gene_exon_dd=1;
      else flag_gene_exon_dd=0;
      end;
   else do;
       if flag_gene_cd4_on=0 and flag_gene_cd8_on=0 and flag_gene_cd19_on=0 then flag_gene_cell_specific=.;
       else flag_gene_cell_specific=1;
   end;
run;

proc freq data=flag_gene_dd;
   tables flag_gene_cell_specific*flag_gene_monoexon flag_gene_exon_dd*flag_gene_monoexon;
run;

/*
flag_gene_cell_specific
          flag_gene_monoexon

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|       1|  Total
---------+--------+--------+
       0 |   9441 |  13795 |  23236
         |  33.73 |  49.28 |  83.01
         |  40.63 |  59.37 |
         |  91.69 |  77.96 |
---------+--------+--------+
       1 |    856 |   3901 |   4757
         |   3.06 |  13.94 |  16.99
         |  17.99 |  82.01 |
         |   8.31 |  22.04 |
---------+--------+--------+
Total       10297    17696    27993
            36.78    63.22   100.00

     Frequency Missing = 21721

           The SAS System         17:39 W

   flag_gene_exon_dd
             flag_gene_monoexon

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 |   5893 |  13795 |  19688
            |  25.36 |  59.37 |  84.73
            |  29.93 |  70.07 |
            |  62.42 | 100.00 |
   ---------+--------+--------+
          1 |   3548 |      0 |   3548
            |  15.27 |   0.00 |  15.27
            | 100.00 |   0.00 |
            |  37.58 |   0.00 |
   ---------+--------+--------+
   Total        9441    13795    23236
               40.63    59.37   100.00

        Frequency Missing = 26478

*/

/* Add genes that are differentially spliced */

data ds_genes;
   set event.t1d_qs_flag_gene_fdr_all;
   keep gene_id flag_cell_by_fus_fdr05 flag_cell_by_fus_fdr10 flag_cell_by_fus_fdr20;
run;

proc sort data=ds_genes;
  by gene_id;
proc sort data=flag_gene_dd;
  by gene_id;
run;


data flag_gene_dd_ds;
  merge flag_gene_dd (in=in1) ds_genes (in=in2);
  by gene_id;
  if in1;
run;

proc freq data=flag_gene_dd_ds;
   where flag_gene_monoexon=0;
   tables flag_gene_exon_dd*flag_cell_by_fus_fdr05
          flag_gene_exon_dd*flag_cell_by_fus_fdr10
          flag_gene_exon_dd*flag_cell_by_fus_fdr20;
run;

proc freq data=flag_gene_dd_ds noprint;
   tables flag_gene_cell_specific*flag_gene_monoexon*flag_gene_exon_dd*flag_cell_by_fus_fdr05
          *flag_cell_by_fus_fdr10*flag_cell_by_fus_fdr20 / out=gene_count;
run;

proc print data=gene_count;
run;


/* Make permenant */

data event.t1d_flag_gene_dd_ds_exons;
  set flag_gene_dd_ds;
run;



/* 
                   flag_     flag_    flag_cell_   flag_cell_   flag_cell_
   flag_gene_      gene_     gene_      by_fus_      by_fus_      by_fus_
 cell_specific   monoexon   exon_dd      fdr05        fdr10        fdr20     COUNT

       .             0         .           .            .            .       21721
       0             0         0           .            .            .        3227

       0             0         0           0            0            0        1092

       0             0         0           0            0            1         152
       0             0         0           0            1            1         131
       0             0         0           1            1            1        1291

       0             0         1           .            .            .        2446
       0             0         1           0            0            0         270
       0             0         1           0            0            1          62
       0             0         1           0            1            1          28
       0             0         1           1            1            1         742

       0             1         0           .            .            .       13795
       1             0         .           .            .            .         856
       1             1         .           .            .            .        3901


21721 genes off
3901+856 genes that are cell specific

13795 non-cell specific monoexon genes

FDR 5%:
1092+152+131 non-cell specific, multi-exon genes that are not DS or DD
1291 non-cell specific, multi-exon genes that are DS but not DD
2446+270+62+28 genes that are DD but not DS
742 genes both DD and DS

FDR 10%:
1092+152 non-cell specific, multi-exon genes that are not DS or DD
1291+131 non-cell specific, multi-exon genes that are DS but not DD
2446+270+62 genes that are DD but not DS
742+28 genes both DD and DS

FDR 20%:
1092 non-cell specific, multi-exon genes that are not DS or DD
1291+152+131 non-cell specific, multi-exon genes that are DS but not DD
2446+270 genes that are DD but not DS
742+62+28  genes both DD and DS

*/


