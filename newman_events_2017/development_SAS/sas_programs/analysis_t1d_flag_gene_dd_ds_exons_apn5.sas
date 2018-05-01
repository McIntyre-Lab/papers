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
   set event.t1d_flag_gene_on_by_cell_apn5;
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

OFF
      .             0         .           .            .            .       35431

CELL SPEC
      1             0         .           .            .            .         102
      1             1         .           .            .            .         195

MONOEXON
      0             1         0           .            .            .        1831

DD
      0             0         1           0            0            0           4
      0             0         1           0            0            1           1
DS
      0             0         0           1            1            1        1727
      0             1         0           1            1            1         252

DS+DD
      0             0         1           1            1            1          54

None
      0             0         0           0            0            0        1179
      0             0         0           0            0            1         185
      0             0         0           0            1            1         136
      0             1         0           0            0            0         179
      0             1         0           0            0            1          28
      0             1         0           0            1            1          23





35431 genes off
102+195=297 genes that are cell specific

1831 non-cell specific monoexon genes

FDR 5%:
1179+185+136+179+28+23=1730 non-cell specific, multi-exon genes that are not DS or DD
1727+252=1979 non-cell specific, multi-exon genes that are DS but not DD
4+1=5 that are DD but not DS
54 genes both DD and DS

1979+5+54=2038 total AS

*/


