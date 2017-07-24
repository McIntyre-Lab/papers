ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Add in RefSeq ID to PacBio transcripts. I want to carry forward only the PB genes that have RefSeq IDs */

data pacbio_as;
  set event.pacbio_flag_gene_alt_spliced;
  drop _TYPE_ _FREQ_;
  num_pb_genes=1;
run;

data pb2gene;
  set event.pacbio_gene_to_refseq;
  where gene_id ne "NovelGene";
  keep pacbio_gene_id gene_id;
run;

proc sort data=pb2gene;
   by pacbio_gene_id;
proc sort data=pacbio_as;
   by pacbio_gene_id;
run;

data pb2gene_as;
  merge pb2gene (in=in1) pacbio_as (in=in2);
  by pacbio_gene_id;
  if in1 and in2;
run;

proc freq data=pb2gene_as noprint;
  tables pacbio_gene_id / out=gene_count;
run;

data pb2drop;
  set gene_count;
  where count > 1;
  keep pacbio_gene_id;
run;

proc sort data=pb2gene_as;
  by pacbio_gene_id;
proc sort data=pb2drop;
  by pacbio_gene_id;
run;

data pb2gene_as2;
  merge pb2gene_as (in=in1) pb2drop (in=in2);
  by pacbio_gene_id;
  if in2 then delete;
run; *6700 genes kept;

proc sort data=pb2gene_as2;
   by gene_id;
proc means data=pb2gene_as2 noprint;
   by gene_id;
   var num_pb_genes num_isoforms_total num_isoforms_both num_isoforms_old  num_isoforms_npc
       flag_pb_single_iso_gene flag_pb_gene_as flag_pb_gene_specific;
   output out=sum_pb_flags sum=;
run;

/* Redo some flags */

data pb2gene_as3;
  set sum_pb_flags;
   if num_isoforms_total=1 then flag_pb2rs_single_iso_gene=1; else flag_pb2rs_single_iso_gene=0;
   if num_isoforms_both > 0 and num_isoforms_npc > 0 and num_isoforms_old > 0 then do;
          flag_pb2rs_gene_as=1;
          flag_pb2rs_gene_specific=0;
          flag_pb2rs_gene_off=0;
          end;
   if num_isoforms_both > 0 and num_isoforms_npc = 0 and num_isoforms_old > 0 then do;
          flag_pb2rs_gene_as=1;
          flag_pb2rs_gene_specific=0;
          flag_pb2rs_gene_off=0;
          end;
   if num_isoforms_both > 0 and num_isoforms_npc > 0 and num_isoforms_old = 0 then do;
          flag_pb2rs_gene_as=1;
          flag_pb2rs_gene_specific=0;
          flag_pb2rs_gene_off=0;
          end;
   if num_isoforms_both > 0 and num_isoforms_npc = 0 and num_isoforms_old = 0 then do;
          flag_pb2rs_gene_as=0;
          flag_pb2rs_gene_specific=0;
          flag_pb2rs_gene_off=0;
          end;
   if num_isoforms_both = 0 and num_isoforms_npc > 0 and num_isoforms_old > 0 then do;
          flag_pb2rs_gene_as=1;
          flag_pb2rs_gene_specific=0;
          flag_pb2rs_gene_off=0;
          end;
   if num_isoforms_both = 0 and num_isoforms_npc = 0 and num_isoforms_old > 0 then do;
          flag_pb2rs_gene_as=0;
          flag_pb2rs_gene_specific=1;
          flag_pb2rs_gene_off=0; end;
   if num_isoforms_both = 0 and num_isoforms_npc > 0 and num_isoforms_old = 0 then do;
          flag_pb2rs_gene_as=0;
          flag_pb2rs_gene_specific=1;
          flag_pb2rs_gene_off=0;
          end;
   if num_isoforms_both = 0 and num_isoforms_npc = 0 and num_isoforms_old = 0 then do;
          flag_pb2rs_gene_as=0;
          flag_pb2rs_gene_specific=0;
          flag_pb2rs_gene_off=1;
          end;
run;

/* Count */

proc freq data=pb2gene_as3 noprint;
  tables flag_pb2rs_single_iso_gene*flag_pb2rs_gene_as*flag_pb2rs_gene_specific*flag_pb2rs_gene_off / out=iso_count;
proc print data=iso_count;
run;

/*
 flag_pb2rs_     flag_
 single_iso_     pb2rs_     flag_pb2rs_     flag_pb2rs_
     gene       gene_as    gene_specific      gene_off     COUNT

      0            0             0               0          2423
      0            0             1               0            47
      0            1             0               0           677
      1            0             0               0          3195
      1            0             1               0           140


Converted to Refseq IDs:
3195 genes with 1 isoform
140 genes with 1 isoform that have cell-specific expression

2423 genes that have multiple isoforms and expressed in both cell types and not AS
47 genes that have multi isoforms and are cell specific
677 genes that have multi isoforms and are AS

677/(677+2423)=22% of multi-PB isoform RefSeq genes in both cell types have AS
Or of all genes, 10%
*/

/* Make permenant */

data event.pb2refseq_flag_gene_alt_spliced;
   set pb2gene_as3;
run;

