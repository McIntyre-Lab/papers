ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* What is the overlap between genes with DD events and genes with DD PB isoforms?

   (1) If a gene only has exonic sequence that a cell-specific then gene is cell-specific (ignore junctions here)
   (2) If a gene has either a fragment OR a junction that is DD, then gene is AS

   Essentially, use fragments to decide if gene is cell-specific or not, and if it has only one exon or not
   then use the combination of fragments and junctions for determining AS

*/

* flag Event Gene if only one annotated transcript;

data gene2xs;
  set event.feature2xs2gene_exp_only_nomulti;
  keep transcript_id gene_id;
run;

proc sort data=gene2xs nodup;
  by gene_id transcript_id;
proc freq data=gene2xs noprint;
  tables gene_id / out =xs_count;
run;

data flag_1iso_gene;
  set xs_count;
  if count=1 then flag_single_xs_gene=1;
  else flag_single_xs_gene=0;
  keep gene_id flag_single_xs_gene;
run;

data event_as;
  set event.genes_all_events_on_by_cell;
  flag_gene_single_exon=flag_single_fusion_gene;
  flag_gene_specific=flag_gene_fusions_specific;
  flag_gene_off=flag_gene_fusions_off;
  if sum(num_fusions_old,num_fusions_npc,num_fusions_both)=1 then flag_gene_single_exon_dtct=1;
  else flag_gene_single_exon_dtct=0;
  if flag_Gene_specific=0 and flag_gene_junc_specific=1 then flag_gene_as=1;
  else if flag_gene_fragments_as=1 or flag_gene_junc_as=1 then flag_gene_as=1;
  else flag_gene_as=0;
  keep gene_id flag_gene_single_exon flag_gene_off flag_gene_specific flag_gene_single_exon_dtct flag_gene_as;
run;

proc sort data=event_as;
  by gene_id;
proc sort data=flag_1iso_gene;
  by gene_id;
run;

data event_as2;
  merge event_as (in=in1) flag_1iso_gene (in=in2);
  by gene_id;
  if in1 and in2;
run;


data pb2refseq_as1;
   set event.pb2refseq_flag_gene_alt_spliced;
   keep gene_id flag_pb2rs_single_iso_gene flag_pb2rs_gene_as flag_pb2rs_gene_specific;
run;


data pb2refseq_as2;
   set event.pb2refseq_flag_gene_alt_spliced;
   if num_pb_genes > 1 then delete;
   keep gene_id flag_pb_single_iso_gene flag_pb_gene_as flag_pb_gene_specific;
run;

proc sort data=event_as2;
  by gene_id;
proc sort data=pb2refseq_as1;
  by gene_id;
proc sort data=pb2refseq_as2;
  by gene_id;
run;

data overlap1;
  merge event_as2 (in=in1) pb2refseq_as1 (in=in2);
  by gene_id;
  if in1 and in2;
run; *4719 genes in common;

data overlap2;
  merge event_as2 (in=in1) pb2refseq_as2 (in=in2);
  by gene_id;
  if in1 and in2;
run;  *4561 genes in common;

proc freq data=overlap1 noprint;
  tables flag_gene_single_exon*flag_gene_off*flag_gene_specific*flag_single_xs_gene*
         flag_gene_single_exon_dtct*flag_pb2rs_single_iso_gene*
         flag_pb2rs_gene_specific*flag_gene_as*flag_pb2rs_gene_as / out=pb2event_check1;
run;

proc export data=pb2event_check1 outfile="!MCLAB/event_analysis/analysis_output/pacbio_vs_events_alt_splicing_overlap_no_pb2multi_refseq.csv"
  dbms=csv replace;
run;


proc freq data=overlap2 noprint;
  tables flag_gene_single_exon*flag_gene_off*flag_gene_specific*flag_single_xs_gene*
         flag_gene_single_exon_dtct*flag_pb_single_iso_gene*
         flag_pb_gene_specific*flag_gene_as*flag_pb_gene_as / out=pb2event_check2;
run;

proc export data=pb2event_check1 outfile="!MCLAB/event_analysis/analysis_output/pacbio_vs_events_alt_splicing_overlap_no_refseq2multi_pb.csv"
  dbms=csv replace;
run;


