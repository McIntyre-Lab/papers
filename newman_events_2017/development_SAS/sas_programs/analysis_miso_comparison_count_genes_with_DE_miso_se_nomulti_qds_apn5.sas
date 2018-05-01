ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Check that the set of genes with DE MISO SE are the same set with DE/DD SE in Event analysis */



/* Get list of genes with MISO events and ES junctions detected/analyzed */

data genes_to_keep;
  set event.miso_refseq_exnskp_cmpr_dtct_qds;
  where flag_has_refseq=1 and flag_has_miso_se_dtct=1 ;
  keep gene_id flag_has_multiple_exons;
run;

/* Merge with MISO */

data miso_diff_se;
   set event.miso_diff_se_refseq;
  where flag_refseq_match=1;
run;


data genes2keep;
   set event.mm10_flag_gene_dd_ds_exons_apn5;
   if flag_gene_cell_specific ne 0 then delete;
   if flag_gene_monoexon=1 and flag_cell_by_fus_fdr05=. then delete;
run;


data dd_genes;
   set genes2keep;
   if flag_gene_exon_dd=1;
   keep gene_id;
run; *1560 genes;


data qds_genes;
   set genes2keep;
   if flag_cell_by_fus_fdr05=1;
   keep gene_id;
run; *153 genes;


proc sort data=miso_diff_se nodup;
   by gene_id;
proc sort data=dd_genes nodup;
   by gene_id;
proc sort data=qds_genes nodup;
   by gene_id;
proc sort data=genes_to_keep nodup;
   by gene_id;
run;

data miso_to_exonskip;
  merge miso_diff_se (in=in1) dd_genes (in=in2) qds_genes (in=in3);
  by gene_id;
  if not in1 then do;
      num_diff_se_bf10=0;
      num_diff_se_bf5=0; end;
  if in2 then flag_gene_dd=1; else flag_gene_dd=0;
  if in3 then flag_gene_qds=1; else flag_gene_qds=0;
  if in1 then output;
run;

data miso_to_exonskip2;
  merge miso_to_exonskip (in=in1) genes_to_keep (in=in2);
  by gene_id;
  if in1 and in2;
run; *654 genes;

data flag_genes;
  set miso_to_exonskip2;
  if num_diff_se_bf10 > 0 then flag_miso_se_diff_bf10=1; else flag_miso_se_diff_bf10=0;
  if num_diff_se_bf5 > 0 then flag_miso_se_diff_bf5=1; else flag_miso_se_diff_bf5=0;
run;

data event.miso_refseq_exonskip_cmpr_dd_qds;
   set flag_genes;
run;


proc freq data=event.miso_refseq_exonskip_cmpr_dd_qds noprint;
   tables flag_miso_se_diff_bf10*flag_gene_dd*flag_gene_qds / out=gene_count_bf10;
   tables flag_miso_se_diff_bf5*flag_gene_dd*flag_gene_qds / out=gene_count_bf5;
run;

proc print data=gene_count_bf10;
run;

proc print data=gene_count_bf5;
run;



/*
BF>10:
 flag_miso_
  se_diff_      flag_       flag_
    bf10       gene_dd    gene_qds    COUNT

      0           0           0        428
      0           0           1          6
      0           1           0        195
      0           1           1         19
      1           0           0          4
      1           1           0          1
      1           1           1          1


BF>5:
flag_miso_
 se_diff_      flag_       flag_
    bf5       gene_dd    gene_qds    COUNT

     0           0           0        426
     0           0           1          6
     0           1           0        194
     0           1           1         18
     1           0           0          6
     1           1           0          2
     1           1           1          2


*/


