ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Comparison of alternatively spliced genes between Event Analysis and PacBio

I want to compare:
(1) Genes that have multiple exons, have exons in both NPCs and OLD and have exons in either NPCs or OLDs
- to -
(2) Genes with multiple PB isoforms, have PB isoforms in both NPCs and OLDs and have isoforms in either NPCs or OLDs

This is the "Event Analysis AS genes" vs. "PacBio AS genes" comparison

Genes that have only one exon cannot be alternatively spliced
Genes that are only detected in one cell type are cell-specific genes and not AS

PB genes with only 1 isoform cannot be AS
PB genes with isoforms only in one cell type are cell-specific and not AS

*/

/* For PacBio genes, I want to count isoforms expressed in:
(1) Both NPCs and OLDs
(2) Only NPCs
(3) Only OLDs

Then, I want to flag:
(1) If gene only has one isoform (thus cannot be AS)
(2) Gene has cell-specific expression (thus cannot be AS)
(3) Gene is alternatively spliced 

*/

data flag_pb_on;
  length pacbio_gene_id $15.;
  set event.pacbio_xscripts_flag_cell_exp;
  pacbio_gene_id=catt("PB.",scan(pacbio_id,2,"."));

  if flag_npc_expressed=1 and flag_old_expressed=1 then flag_isoform_both=1; else flag_isoform_both=0;
  if flag_npc_expressed=1 and flag_old_expressed=0 then flag_isoform_npc=1; else flag_isoform_npc=0;
  if flag_npc_expressed=0 and flag_old_expressed=1 then flag_isoform_old=1; else flag_isoform_old=0;
  flag_isoform=1;
run;

proc sort data=flag_pb_on;
  by pacbio_gene_id;
proc means data=flag_pb_on noprint;
  by pacbio_gene_id;
  var flag_isoform flag_isoform_both flag_isoform_npc flag_isoform_old;
  output out=pb_iso_count sum(flag_isoform)=num_isoforms_total
                          sum(flag_isoform_both)=num_isoforms_both
                          sum(flag_isoform_npc)=num_isoforms_npc
                          sum(flag_isoform_old)=num_isoforms_old;
run;

/* Flag if PacBio genes are alternatively spliced, have only 1 isoform, and have cell-specific expression */

data flag_pb_as;
   set pb_iso_count;
   if num_isoforms_total=1 then flag_pb_single_iso_gene=1; else flag_pb_single_iso_gene=0;
   if num_isoforms_both > 0 and num_isoforms_npc > 0 and num_isoforms_old > 0 then do;
          flag_pb_gene_as=1;
          flag_pb_gene_specific=0;
          flag_pb_gene_off=0;
          end;
   if num_isoforms_both > 0 and num_isoforms_npc = 0 and num_isoforms_old > 0 then do;
          flag_pb_gene_as=1;
          flag_pb_gene_specific=0;
          flag_pb_gene_off=0;
          end;
   if num_isoforms_both > 0 and num_isoforms_npc > 0 and num_isoforms_old = 0 then do;
          flag_pb_gene_as=1;
          flag_pb_gene_specific=0;
          flag_pb_gene_off=0;
          end;
   if num_isoforms_both > 0 and num_isoforms_npc = 0 and num_isoforms_old = 0 then do;
          flag_pb_gene_as=0;
          flag_pb_gene_specific=0;
          flag_pb_gene_off=0;
          end;
   if num_isoforms_both = 0 and num_isoforms_npc > 0 and num_isoforms_old > 0 then do;
          flag_pb_gene_as=1;
          flag_pb_gene_specific=0;
          flag_pb_gene_off=0;
          end;
   if num_isoforms_both = 0 and num_isoforms_npc = 0 and num_isoforms_old > 0 then do;
          flag_pb_gene_as=0;
          flag_pb_gene_specific=1;
          flag_pb_gene_off=0; end;
   if num_isoforms_both = 0 and num_isoforms_npc > 0 and num_isoforms_old = 0 then do;
          flag_pb_gene_as=0;
          flag_pb_gene_specific=1;
          flag_pb_gene_off=0;
          end;
   if num_isoforms_both = 0 and num_isoforms_npc = 0 and num_isoforms_old = 0 then do;
          flag_pb_gene_as=0;
          flag_pb_gene_specific=0;
          flag_pb_gene_off=1;
          end;
run;

/* Count */

proc freq data=flag_pb_as noprint;
  tables flag_pb_single_iso_gene*flag_pb_gene_as*flag_pb_gene_specific*flag_pb_gene_off / out=iso_count;
proc print data=iso_count;
run;


/*
 flag_pb_
  single_    flag_pb_    flag_pb_gene_    flag_pb_
 iso_gene     gene_as       specific      gene_off    COUNT    PERCENT

     0           0             0              0        2526    32.8949
     0           0             1              0          53     0.6902
     0           1             0              0         700     9.1158
     1           0             0              0        4152    54.0695
     1           0             1              0         248     3.2296

248 genes have a single isoform that is only detected in one cell type
4152 genes have a single isoform that is detected in both cell types

2526 genes that have multiple isoforms, but do not have cell-specific expression or AS
53 genes have multiple isoforms that are all detected in one cell type
700 genes that have multiple isoforms and have evidence of alternative splicing

Thus 700/(700+2526) = 22% of genes with multiple isoforms and expressed in both cell types have AS
Or if all genes then 9%

*/


/* Make permenant -- next step will be to merge in RefSeq IDs */

data event.pacbio_flag_gene_alt_spliced;
  set flag_pb_as;
run;



