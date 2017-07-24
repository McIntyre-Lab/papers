ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Process PacBio data: we only want to look at genes and transcripts with a RefSeq match.
   We also want to remove genes that aren't expressed, or have multigene components
   such that the comparison between Event Analysis and PacBio is as "clean" as possible

   For now, I am going to remove genes with multigene components and leave the "gene expressed"
   flag in the PB2RefSeq dataset */

/* Get list of expressed genes */

data genes_exp;
   set event.flag_gene_expressed;
   keep gene_id;
run;

/* PacBio to Refseq gene */

data pb2refseq;
   set event.pacbio2refseq_gene;
run;

proc sort data=pb2refseq;
   by gene_id;
proc sort data=genes_exp;
   by gene_id;
run;

data pb2refseq_filtered no_pb no_rs;
   merge genes_exp (in=in1) pb2refseq (in=in2);
   by gene_id;
   if in1 and in2 then output pb2refseq_filtered;
   else if in1 then output no_pb;
   else output no_rs;
run;

*6926 PB transcripts with RefSeq IDs kept;
*23698 transcripts without multigene pieces with no matching PB gene;
*2727 PB transcripts without a RefSeq ID;

/* PacBio to Refseq transcripts */

data pb2xs;
   set event.pacbio2refseq_id;
run;

data pb_xs_kept;
   set pb2refseq_filtered;
   keep pacbio_id;
run;

proc sort data=pb2xs;
   by pacbio_id;
proc sort data=pb_xs_kept nodup;
   by pacbio_id;
run;

data pb2xs_filtered;
  merge pb2xs (in=in1) pb_xs_kept (in=in2);
  if in1 and in2;
run;

/* All PacBio transcripts */

data pb_gene2keep;
   set pb2refseq_filtered;
   keep pacbio_gene_id;
run;

data pb_xs;
   set event.pacbio_transcripts;
run;

proc sort data=pb_gene2keep nodup; 
  by pacbio_gene_id;
proc sort data=pb_xs;
  by pacbio_gene_id;
run;

data pb_xs_filtered;
  merge pb_gene2keep (in=in1) pb_xs (in=in2);
  by pacbio_gene_id;
  if in1 and in2;
run;

/* Make all permenant */

data event.pacbio2refseq_gene_nomulti;
   set pb2refseq_filtered;
run;

data event.pacbio2refseq_id_nomulti;
   set pb2xs_filtered;
run;

data event.pacbio_transcripts_nomulti;
   set pb_xs_filtered;
run;

