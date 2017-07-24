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
   set event.mm10_flag_gene_on_by_cell_apn5;
   /* If there are more exons detected in either cell type than in both cell types together,
      then the gene has differentially detected exons */
    if sum((flag_fusion_nsc_on-flag_fusion_both_on),
             (flag_fusion_old_on-flag_fusion_both_on),
             flag_fusion_both_on) = 1 then flag_gene_monoexon=1;
    else flag_gene_monoexon=0;

   if flag_gene_both_on=1 then do;
      flag_gene_cell_specific=0;
      if flag_fusion_nsc_on > flag_fusion_both_on or flag_fusion_old_on > flag_fusion_both_on then flag_gene_exon_dd=1;
       else flag_gene_exon_dd=0;
      end;
   else do;
       if flag_gene_nsc_on=0 and flag_gene_old_on=0 then flag_gene_cell_specific=.;
       else flag_gene_cell_specific=1;
   end;
run;

proc freq data=flag_gene_dd;
   tables flag_gene_cell_specific*flag_gene_monoexon flag_gene_exon_dd*flag_gene_monoexon;
run;

/*
               flag_gene_monoexon

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |   5724 |    940 |   6664
              |  69.71 |  11.45 |  81.16
              |  85.89 |  14.11 |
              |  86.24 |  59.72 |
     ---------+--------+--------+
            1 |    913 |    634 |   1547
              |  11.12 |   7.72 |  18.84
              |  59.02 |  40.98 |
              |  13.76 |  40.28 |
     ---------+--------+--------+
     Total        6637     1574     8211
                 80.83    19.17   100.00

          Frequency Missing = 20106

     flag_gene_exon_dd
               flag_gene_monoexon

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |   4164 |    940 |   5104
              |  62.48 |  14.11 |  76.59
              |  81.58 |  18.42 |
              |  72.75 | 100.00 |
     ---------+--------+--------+
            1 |   1560 |      0 |   1560
              |  23.41 |   0.00 |  23.41
              | 100.00 |   0.00 |
              |  27.25 |   0.00 |
     ---------+--------+--------+
     Total        5724      940     6664
                 85.89    14.11   100.00

          Frequency Missing = 21653



*/

/* Add genes that are differentially spliced */

data ds_genes;
   set event.mm10_NPCvOLD_quantspl_flag_fdr;
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

data event.mm10_flag_gene_dd_ds_exons_apn5;
  set flag_gene_dd_ds;
run;



/* 
                  flag_     flag_    flag_cell_   flag_cell_   flag_cell_
  flag_gene_      gene_     gene_      by_fus_      by_fus_      by_fus_
cell_specific   monoexon   exon_dd      fdr05        fdr10        fdr20     COUNT   PERCENT

      .             0         .           .            .            .       20106     .

      0             0         0           0            0            0        4059   64.7162
      0             0         0           0            0            1          30    0.4783
      0             0         0           0            1            1          21    0.3348

      0             0         0           1            1            1          54    0.8610

      0             0         1           0            0            0        1388   22.1301
      0             0         1           0            0            1          51    0.8131
      0             0         1           0            1            1          28    0.4464

      0             0         1           1            1            1          93    1.4828

      0             1         0           0            0            0         527    8.4024
      0             1         0           0            0            1          10    0.1594
      0             1         0           0            1            1           5    0.0797

      0             1         0           1            1            1           6    0.0957

      0             1         0           .            .            .         392     .


      1             0         .           .            .            .         913     .
      1             1         .           .            .            .         634     .


20106 genes off
913+634=1547 genes that are cell specific
392 non-cell specific monoexon genes

FDR 5%:
4059+30+21+527+10+5=4652 non-cell specific, multi-exon genes that are not DS or DD

54+6=60 non-cell specific, multi-exon genes that are DS but not DD

1388+51+28=1560 genes that are DD but not DS

93 genes both DD and DS


*/


