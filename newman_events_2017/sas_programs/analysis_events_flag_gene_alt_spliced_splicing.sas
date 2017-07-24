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

/* For each gene, count junctions that are on in both, NPCs only, OLDs only, or off */

data flag_event_dd;
   set event.flag_splicing_on;
   if flag_event_nsc_on=1 and flag_event_old_on=1 then flag_event_both=1; else flag_event_both=0;
   if flag_event_nsc_on=1 and flag_event_old_on=0 then flag_event_npc=1; else flag_event_npc=0;
   if flag_event_nsc_on=0 and flag_event_old_on=1 then flag_event_old=1; else flag_event_old=0;
   if flag_event_nsc_on=0 and flag_event_old_on=0 then flag_event_off=1; else flag_event_off=0;
   keep event_id flag_event_both flag_event_npc flag_event_old flag_event_off;
run;

data event2keep;
   set event.flagged_event_length;
   where flag_event_short=0;
   keep event_id;
run;

proc sort data=event2keep;
  by event_id;
proc sort data=flag_event_dd;
  by event_id;
run;

data flag_event_dd2;
   merge flag_event_dd (in=in1) event2keep (in=in2);
   by event_id;
   if in1 and in2;
run;


data event2gene;
  set evspl.splicing_events_annot_refseq;
  keep event_id gene_id;
run;


data gene_exp;
   set event.flag_gene_expressed;
   where flag_gene_expressed=1;
   keep gene_id;
run;

proc sort data=event2gene nodup;
   by gene_id event_id;
proc sort data=gene_exp;
   by gene_id;
run;

data event2gene2;
  merge event2gene (in=in1) gene_exp (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=event2gene2;
   by event_id;
proc sort data=flag_event_dd2;
   by event_id;
run;

data flag_event_dd3;
  merge flag_event_dd2 (in=in1) event2gene2 (in=in2);
  by event_id;
  if in1 and in2;
  flag_event=1;
run;

proc sort data=flag_event_dd3;
   by gene_id;
proc means data=flag_event_dd3 noprint;
   by gene_id;
   var flag_event_both flag_event_npc flag_event_old flag_event_off flag_event;
   output out=junc_on_counts sum(flag_event_both)=num_junc_both
                             sum(flag_event_npc)=num_junc_npc
                             sum(flag_event_old)=num_junc_old
                             sum(flag_event_off)=num_junc_off
                             sum(flag_event)=num_junc_total;
run;


data check;
  set junc_on_counts;
  if num_junc_both+num_junc_npc+num_junc_old+num_junc_off = num_junc_total then delete;
run;

/* Make permenant */

data event.genes_junctions_on_by_cell;
  set junc_on_counts;
  drop _TYPE_ _FREQ_;
run;


