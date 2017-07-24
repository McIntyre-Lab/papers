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

/* Merge fragment/fusion/splicing counts and flag if:
-gene is cell-specific
-gene has only 1 exon
-gene is alternatively spliced
*/

data frags;
   set event.genes_frags_on_by_cell;
run;
data fusions;
   set event.genes_fusions_on_by_cell;
run;
data juncs;
   set event.genes_junctions_on_by_cell;
run;

proc sort data=frags;
  by gene_id;
proc sort data=fusions;
  by gene_id;
proc sort data=juncs;
  by gene_id;
run;

data all_counts_by_gene;
  merge  fusions frags juncs ;
   by gene_id;
run;

/* Flag if gene is cell-specific, 1 exon, is AS, etc. */

data flag_gene_exp;
   set all_counts_by_gene;
   /* Flag if single exon gene */
   if num_fusions_total=1 then flag_single_fusion_gene=1; else flag_single_fusion_gene=0;
   if num_fragments_total=1 then flag_single_fragment_gene=1; else flag_single_fragment_gene=0;
run;

data flag_gene_dd_fusions;
   set flag_gene_exp;
  /* Flag if gene expression is cell specific or if alternatively spliced */
  *fusions;
  if num_fusions_both > 0 and num_fusions_npc > 0 and num_fusions_old > 0 then do;
       flag_gene_fusions_as=1;
       flag_gene_fusions_specific=0;
       flag_gene_fusions_off=0;
       end;
  if num_fusions_both > 0 and num_fusions_npc = 0 and num_fusions_old > 0 then do;
       flag_gene_fusions_as=1;
       flag_gene_fusions_specific=0;
       flag_gene_fusions_off=0;
       end;
  if num_fusions_both > 0 and num_fusions_npc > 0 and num_fusions_old = 0 then do;
       flag_gene_fusions_as=1;
       flag_gene_fusions_specific=0;
       flag_gene_fusions_off=0;
       end;
  if num_fusions_both > 0 and num_fusions_npc = 0 and num_fusions_old = 0 then do;
       flag_gene_fusions_as=0;
       flag_gene_fusions_specific=0;
       flag_gene_fusions_off=0;
       end;
  if num_fusions_both = 0 and num_fusions_npc > 0 and num_fusions_old > 0 then do;
       flag_gene_fusions_as=1;
       flag_gene_fusions_specific=0;
       flag_gene_fusions_off=0;
       end;
  if num_fusions_both = 0 and num_fusions_npc = 0 and num_fusions_old > 0 then do;
       flag_gene_fusions_as=0;
       flag_gene_fusions_specific=1;
       flag_gene_fusions_off=0;
       end;
  if num_fusions_both = 0 and num_fusions_npc > 0 and num_fusions_old = 0 then do;
       flag_gene_fusions_as=0;
       flag_gene_fusions_specific=1;
       flag_gene_fusions_off=0;
       end;
  if num_fusions_both = 0 and num_fusions_npc = 0 and num_fusions_old = 0 then do;
       flag_gene_fusions_as=0;
       flag_gene_fusions_specific=0;
       flag_gene_fusions_off=1;
       end;
run;



data flag_gene_dd_frags;
   set flag_gene_dd_fusions;
  /* Flag if gene expression is cell specific or if alternatively spliced */
  *fragments;
  if num_fragments_both > 0 and num_fragments_npc > 0 and num_fragments_old > 0 then do;
       flag_gene_fragments_as=1;
       flag_gene_fragments_specific=0;
       flag_gene_fragments_off=0;
       end;
  if num_fragments_both > 0 and num_fragments_npc = 0 and num_fragments_old > 0 then do;
       flag_gene_fragments_as=1;
       flag_gene_fragments_specific=0;
       flag_gene_fragments_off=0;
       end;
  if num_fragments_both > 0 and num_fragments_npc > 0 and num_fragments_old = 0 then do;
       flag_gene_fragments_as=1;
       flag_gene_fragments_specific=0;
       flag_gene_fragments_off=0;
       end;
  if num_fragments_both > 0 and num_fragments_npc = 0 and num_fragments_old = 0 then do;
       flag_gene_fragments_as=0;
       flag_gene_fragments_specific=0;
       flag_gene_fragments_off=0;
       end;
  if num_fragments_both = 0 and num_fragments_npc > 0 and num_fragments_old > 0 then do;
       flag_gene_fragments_as=1;
       flag_gene_fragments_specific=0;
       flag_gene_fragments_off=0;
       end;
  if num_fragments_both = 0 and num_fragments_npc = 0 and num_fragments_old > 0 then do;
       flag_gene_fragments_as=0;
       flag_gene_fragments_specific=1;
       flag_gene_fragments_off=0;
       end;
  if num_fragments_both = 0 and num_fragments_npc > 0 and num_fragments_old = 0 then do;
       flag_gene_fragments_as=0;
       flag_gene_fragments_specific=1;
       flag_gene_fragments_off=0;
       end;
  if num_fragments_both = 0 and num_fragments_npc = 0 and num_fragments_old = 0 then do;
       flag_gene_fragments_as=0;
       flag_gene_fragments_specific=0;
       flag_gene_fragments_off=1;
       end;
run;


data flag_gene_dd_juncs;
   set flag_gene_dd_frags;
  /* Flag if gene expression is cell specific or if alternatively spliced */
  *splicing ;
   if num_junc_total > 0 then do;
  if num_junc_both > 0 and num_fragments_npc > 0 and num_fragments_old > 0 then do;
       flag_gene_junc_as=1;
       flag_gene_junc_specific=0;
       flag_gene_junc_off=0;
       end;
  if num_junc_both > 0 and num_junc_npc = 0 and num_junc_old > 0 then do;
       flag_gene_junc_as=1;
       flag_gene_junc_specific=0;
       flag_gene_junc_off=0;
       end;
  if num_junc_both > 0 and num_junc_npc > 0 and num_junc_old = 0 then do;
       flag_gene_junc_as=1;
       flag_gene_junc_specific=0;
       flag_gene_junc_off=0;
       end;
  if num_junc_both > 0 and num_junc_npc = 0 and num_junc_old = 0 then do;
       flag_gene_junc_as=0;
       flag_gene_junc_specific=0;
       flag_gene_junc_off=0;
       end;
  if num_junc_both = 0 and num_junc_npc > 0 and num_junc_old > 0 then do;
       flag_gene_junc_as=1;
       flag_gene_junc_specific=0;
       flag_gene_junc_off=0;
       end;
  if num_junc_both = 0 and num_junc_npc = 0 and num_junc_old > 0 then do;
       flag_gene_junc_as=0;
       flag_gene_junc_specific=1;
       flag_gene_junc_off=0;
       end;
  if num_junc_both = 0 and num_junc_npc > 0 and num_junc_old = 0 then do;
       flag_gene_junc_as=0;
       flag_gene_junc_specific=1;
       flag_gene_junc_off=0;
       end;
  if num_junc_both = 0 and num_junc_npc = 0 and num_junc_old = 0 then do;
       flag_gene_junc_as=0;
       flag_gene_junc_specific=0;
       flag_gene_junc_off=1;
       end;
  end;
run;


/* Count genes that are cell-specifically expressed or alternatively spliced */

proc freq data=flag_gene_dd_juncs noprint;
   tables flag_single_fusion_gene*flag_gene_fusions_as*flag_gene_fusions_specific*flag_gene_fusions_off / out=gene_counts_fus;
proc print data=gene_counts_fus;
run;

/*
                               flag_gene_    flag_gene_
 flag_single_    flag_gene_     fusions_      fusions_
  fusion_gene    fusions_as     specific         off       COUNT    PERCENT

       0              0             0             0         8600    41.1976
       0              0             1             0         2469    11.8275
       0              1             0             0         7975    38.2036
       1              0             0             0         1429     6.8455
       1              0             1             0          402     1.9257

402 genes are single-exon genes and have cell-specific expression
1429 genes are single-exon genes but are expressed in both cell type
2469 genes have multiple exons and have have cell-specific expression

8600 genes are expressed in both cell types but are not alternatively spliced
7975 genes are expressed in both cell types and are alternatively spliced

So, based on fusion-level info 7975/(7975+8600) genes (48%) are alternatively spliced
Or if counting all genes, about 38% 

*/


proc freq data=flag_gene_dd_juncs noprint;
   tables flag_single_fragment_gene*flag_gene_fragments_as*flag_gene_fragments_specific*flag_gene_fragments_off / out=gene_counts_frags;
proc print data=gene_counts_frags;
run;

/*

 flag_single_    flag_gene_    flag_gene_    flag_gene_
   fragment_     fragments_    fragments_    fragments_
     gene            as         specific         off       COUNT    PERCENT

       0              0             0             0         7488    35.8707
       0              0             1             0         2483    11.8946
       0              1             0             0         9194    44.0431
       1              0             0             0         1322     6.3329
       1              0             1             0          388     1.8587

388 genes are single-fragment genes and have cell-specific expression
1322 genes are single-fragment genes but are expressed in both cell type
2483 genes have multiple fragments and  have cell-specific expression

7488 genes are expressed in both cell types but are not alternatively spliced
9194 genes are expressed in both cell types and are alternatively spliced

So, based on fragment-level info 9194/(9194+7488) genes (55%) are alternatively spliced
Or if counting all genes, about 44% of genes
*/


proc freq data=flag_gene_dd_juncs noprint;
   tables flag_single_fragment_gene*flag_gene_junc_as*flag_gene_junc_specific*flag_gene_junc_off / out=gene_counts_junc;
proc print data=gene_counts_junc;
run;

/*

flag_single_     flag_                        flag_
  fragment_      gene_       flag_gene_       gene_
    gene        junc_as    junc_specific    junc_off    COUNT    PERCENT

      0            .             .              .        4619      .
      0            0             0              0        1413     9.7140
      0            0             0              1        5274    36.2574
      0            0             1              0        2670    18.3556
      0            1             0              0        5189    35.6730
      1            .             .              .        1710      .


4619+1710 genes with no junctions

Of genes with junctions:
2670 have cell-specific expression of junctions
1413 genes are expressed in both cell types but are not alternatively spliced
5274 genes have none of their junctions detected
5189 genes are expressed in both cell types and are alternatively spliced

So, based on fragment-level info 5189/(5189+1413) genes (79%) are alternatively spliced
Or if counting all genes, about 36% of genes
*/

/* Make permenant -- think I will use the fusion-level data only for determining if a gene is AS or not */

data event.genes_all_events_on_by_cell;
   set flag_gene_dd_juncs;
run;



