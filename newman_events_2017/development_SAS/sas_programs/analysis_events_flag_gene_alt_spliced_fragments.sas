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

/* For each gene, count fragments that are on in both, NPCs only, OLDs only, or off */

data flag_frag_dd;
   set event.flag_fragment_on;
   if flag_fragment_nsc_on=1 and flag_fragment_old_on=1 then flag_frag_both=1; else flag_frag_both=0;
   if flag_fragment_nsc_on=1 and flag_fragment_old_on=0 then flag_frag_npc=1; else flag_frag_npc=0;
   if flag_fragment_nsc_on=0 and flag_fragment_old_on=1 then flag_frag_old=1; else flag_frag_old=0;
   if flag_fragment_nsc_on=0 and flag_fragment_old_on=0 then flag_frag_off=1; else flag_frag_off=0;
   keep fragment_id flag_frag_both flag_frag_npc flag_frag_old flag_frag_off;
run;

data frag2keep;
   set event.flagged_fragment_length;
   where flag_fragment_lt_min_bp=0;
   keep fragment_id;
run;

proc sort data=frag2keep;
  by fragment_id;
proc sort data=flag_frag_dd;
  by fragment_id;
run;

data flag_frag_dd2;
   merge flag_frag_dd (in=in1) frag2keep (in=in2);
   by fragment_id;
   if in1 and in2;
run;


data frag2gene;
  set mm10.mm10_fragment2exon2gene;
  keep fragment_id gene_id;
run;

data gene_exp;
   set event.flag_gene_expressed;
   where flag_gene_expressed=1;
   keep gene_id;
run;

proc sort data=frag2gene nodup;
   by gene_id fragment_id;
proc sort data=gene_exp;
   by gene_id;
run;

data frag2gene2;
  merge frag2gene (in=in1) gene_exp (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=frag2gene2;
   by fragment_id;
proc sort data=flag_frag_dd2;
   by fragment_id;
run;

data flag_frag_dd3;
  merge flag_frag_dd2 (in=in1) frag2gene2 (in=in2);
  by fragment_id;
  if in1 and in2;
  flag_frag=1;
run;

proc sort data=flag_frag_dd3;
   by gene_id;
proc means data=flag_frag_dd3 noprint;
   by gene_id;
   var flag_frag_both flag_frag_npc flag_frag_old flag_frag_off flag_frag;
   output out=frag_on_counts sum(flag_frag_both)=num_fragments_both
                             sum(flag_frag_npc)=num_fragments_npc
                             sum(flag_frag_old)=num_fragments_old
                             sum(flag_frag_off)=num_fragments_off
                             sum(flag_frag)=num_fragments_total;
run;

/* Make permenant */

data event.genes_frags_on_by_cell;
  set frag_on_counts;
  drop _TYPE_ _FREQ_;
run;


