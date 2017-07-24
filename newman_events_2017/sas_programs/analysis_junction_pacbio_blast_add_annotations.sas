ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

data pacbio_hits;
   set event.blast_junctions_to_pacbio;
run;


/* Add in annotations -- need flags for annot junc and IR, event size, gene_id */

data event_annot;
  set evspl.splicing_events_annot_refseq;
  keep event_id gene_id;
run;

proc sort data=event_annot;
  by event_id;
proc sort data=pacbio_junc;
  by event_id;
run;

data pacbio_blast_w_annot;
   merge pacbio_junc (in=in1) event_annot (in=in2);
   by event_id;
   if in1 and in2;
run;

data pacbio_blast_w_annot2;
   set pacbio_blast_w_annot;
   length pacbio_gene_id $10.;
   pacbio_gene_id=catt("PB.",scan(pacbio_id,2,"."));
run;

data pb2rs_gene;
   set event.pacbio_gene_to_refseq;
   if gene_id="NovelGene" then delete;
   keep gene_id pacbio_gene_id;
run;

proc sort data=pb2rs_gene nodup;
  by pacbio_gene_id;
proc freq data=pb2rs_gene noprint;
  tables pacbio_gene_id / out=gene_count;
proc sort data=gene_count;
  by descending count;
proc print data=gene_count (obs=1);
run; *15 RefSeq genes per PacBio gene max;

data cat_gene;
  array gene[15] $15.;
  retain gene1-gene15;
  set pb2rs_gene;
  by pacbio_gene_id;
  if first.pacbio_gene_id then do;
   call missing(of gene1-gene15);
   records=0;
  end;
  records + 1;
  gene[records] = gene_id;
  if last.pacbio_gene_id then output;
run;

data cat_gene2;
  set cat_gene;
  length refseq_gene_id $200.;
  refseq_gene_id=catx("|", OF gene1-gene15);
  keep pacbio_gene_id refseq_gene_id;
run;


proc sort data=pacbio_blast_w_annot2;
  by pacbio_gene_id;
run;

data pacbio_blast_w_annot3;
  merge pacbio_blast_w_annot2 (in=in1) cat_gene2 (in=in2);
  by pacbio_gene_id;
  if in1 and in2 then do;
        flag_pb_gene_has_refseq=1;
        output; end;
  else if in1 then do;
        flag_pb_gene_has_refseq=0;
        output; end;
run;

/* Make permenant */

data event.blast_junc_pb_hits_with_annot;
   set pacbio_blast_w_annot3;
run;

