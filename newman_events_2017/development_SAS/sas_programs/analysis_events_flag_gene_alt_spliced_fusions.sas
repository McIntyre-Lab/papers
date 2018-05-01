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

/* Count the number of exonic regions, exon fragments and junctions per gene:
(1) In total
(2) Detected in both
(3) Detected in NPCs only
(4) Detected in OLDs only
*/

/* For each gene, count fusions that are on in both, NPCs only, OLDs only, or off */

data flag_fusion_dd;
   set event.flag_fusion_on;
   if flag_fusion_nsc_on=1 and flag_fusion_old_on=1 then flag_fusion_both=1; else flag_fusion_both=0;
   if flag_fusion_nsc_on=1 and flag_fusion_old_on=0 then flag_fusion_npc=1; else flag_fusion_npc=0;
   if flag_fusion_nsc_on=0 and flag_fusion_old_on=1 then flag_fusion_old=1; else flag_fusion_old=0;
   if flag_fusion_nsc_on=0 and flag_fusion_old_on=0 then flag_fusion_off=1; else flag_fusion_off=0;
   keep fusion_id flag_fusion_both flag_fusion_npc flag_fusion_old flag_fusion_off;
run;

data fusion2keep;
   set event.flagged_fusion_length;
   where flag_fusion_lt_min_bp=0;
   keep fusion_id;
run;

proc sort data=fusion2keep;
  by fusion_id;
proc sort data=flag_fusion_dd;
  by fusion_id;
run;

data flag_fusion_dd2;
   merge flag_fusion_dd (in=in1) fusion2keep (in=in2);
   by fusion_id;
   if in1 and in2;
run;


data fusion2gene;
  set mm10.mm10_refseq_fusion_si_info_v2;
  keep fusion_id primary_gene_id;
  rename primary_gene_id=gene_id;
run;

data gene_exp;
   set event.flag_gene_expressed;
   where flag_gene_expressed=1;
   keep gene_id;
run;

proc sort data=fusion2gene nodup;
   by gene_id fusion_id;
proc sort data=gene_exp;
   by gene_id;
run;

data fusion2gene2;
  merge fusion2gene (in=in1) gene_exp (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=fusion2gene2;
   by fusion_id;
proc sort data=flag_fusion_dd2;
   by fusion_id;
run;

data flag_fusion_dd3;
  merge flag_fusion_dd2 (in=in1) fusion2gene2 (in=in2);
  by fusion_id;
  if in1 and in2;
  flag_fusion=1;
run;

proc sort data=flag_fusion_dd3;
   by gene_id;
proc means data=flag_fusion_dd3 noprint;
   by gene_id;
   var flag_fusion_both flag_fusion_npc flag_fusion_old flag_fusion_off flag_fusion;
   output out=fusion_on_counts sum(flag_fusion_both)=num_fusions_both
                             sum(flag_fusion_npc)=num_fusions_npc
                             sum(flag_fusion_old)=num_fusions_old
                             sum(flag_fusion_off)=num_fusions_off
                             sum(flag_fusion)=num_fusions_total;
run;


data check;
  set fusion_on_counts;
  if num_fusions_both+num_fusions_npc+num_fusions_old+num_fusions_off = num_fusions_total then delete;
run;

/* Make permenant */

data event.genes_fusions_on_by_cell;
  set fusion_on_counts;
  drop _TYPE_ _FREQ_;
run;


