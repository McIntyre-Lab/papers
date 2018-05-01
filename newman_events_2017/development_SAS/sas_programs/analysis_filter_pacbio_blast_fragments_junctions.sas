ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Pre-process PacBio BLAST results:
Remove features not expressed
Remove redundant features
Remove features from genes with multigene components
Remove hits from genes not expressed
Check feature is (at least) going to the correct gene, and flag if not
*/

* features expressed;
data feat_on;
   set event.features_w_annotations_nomulti;
   where flag_feature_on=1;
   keep feature_id;
run;

* NR features;
data feat_nr;
  set event.feature2xs2gene_exp_only_nomulti;
  keep feature_id gene_id;
run;

* Genes to keep;
data genes2keep;
  set event.flag_gene_expressed;
  where flag_gene_expressed=1;
  keep gene_id;
run;

/* PacBio BLAST results */

data pacbio_blast;
   set event.pacbio_blast_results;
run;

proc sort data=feat_on nodup;
   by feature_id;
proc sort data=feat_nr nodup;
   by feature_id gene_id;
proc sort data=genes2keep nodup;
   by gene_id;
proc sort data=pacbio_blast nodup;
   by feature_id;
run;

data pacbio_blast_dtct;
  merge pacbio_blast (in=in1) feat_on (in=in2);
   by feature_id;
   if in1 and in2;
run;

data pacbio_blast_nr;
  merge pacbio_blast_dtct (in=in1) feat_nr (in=in2);
  by feature_id;
  if in1 and in2;
run;

proc sort data=pacbio_blast_nr;
   by gene_id;
run;

data pacbio_blast_w_gene;
  merge pacbio_blast_nr (in=in1) genes2keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* For each feature, I want to add in the list of transcripts assigned to it */

data feat2xs;
   set event.feature2xs2gene_exp_only_nomulti;
   keep feature_id transcript_id;
run;

proc sort data=feat2xs nodup;
   by feature_id transcript_id;
proc freq data=feat2xs noprint;
   tables feature_id / out=xs_per_feat;
proc sort data=xs_per_feat;
   by descending COUNT;
proc print data=xs_per_feat (obs=1);
run; *68 transcripts;

data cat_xscript;
  array xs[68] $ 16.;
  retain xs1-xs68;
  set feat2xs;
  by feature_id;
  if first.feature_id then do;
     call missing(of xs1-xs68);
     records = 0;
  end;
  records + 1;
  xs[records]=transcript_id;
  if last.feature_id then output;
run;

data cat_xscript2;
  set cat_xscript;
  length transcript_list $ 1000.;
         transcript_list= catx("|", OF xs1-xs68);
  keep feature_id  transcript_list;
  run;

proc sort data=cat_xscript2;
   by feature_id;
proc sort data=pacbio_blast_w_gene;
   by feature_id;
run;

data pb_blast_w_xs;
  merge pacbio_blast_w_gene (in=in1) cat_xscript2;
  by feature_id;
  if in1;
run;

/* Now get RefSeq IDs for PacBio transcripts, where applicable */

data pb2gene;
  set event.pacbio2refseq_gene_nomulti;
  keep gene_id pacbio_gene_id;
run;

proc sort data=pb2gene nodup;
   by pacbio_gene_id;
run;

proc sort data=pb2gene nodup;
   by pacbio_gene_id gene_id;
proc freq data=pb2gene noprint;
   tables pacbio_gene_id / out=rs_per_pb;
proc sort data=rs_per_pb;
   by descending COUNT;
proc print data=rs_per_pb (obs=1);
run; *2 genes;

data cat_refseq;
  array rs[2] $ 16.;
  retain rs1-rs2;
  set pb2gene;
  by pacbio_gene_id;
  if first.pacbio_gene_id then do;
     call missing(of rs1-rs2);
     records = 0;
  end;
  records + 1;
  rs[records]=gene_id;
  if last.pacbio_gene_id then output;
run;

data cat_refseq2;
  set cat_refseq;
  length refseq_gene_id $ 50.;
         refseq_gene_id= catx("|", OF rs1-rs2);
  keep pacbio_gene_id  refseq_gene_id;
  run;

data pb_xs;
   set event.pacbio_transcripts_nomulti;
   keep pacbio_gene_id pacbio_id;
run;

proc sort data=cat_refseq2 nodup;
   by pacbio_gene_id;
proc sort data=pb_xs nodup;
   by pacbio_gene_id;
run;

data pb2refseq;
  merge cat_refseq2 (in=in1) pb_xs (in=in2);
  by pacbio_gene_id;
  if in1 and in2;
run;

/* Merge in to PacBio BLAST results */

proc sort data=pb2refseq;
   by pacbio_id;
run;

proc sort data=pb_blast_w_xs;
   by pacbio_id;
run;

data blast_w_pb2refseq;
  merge pb_blast_w_xs (in=in1) pb2refseq (in=in2);
  by pacbio_id;
  if in1 and in2;
run;

/* Now we have the set of features/events that are detected, are from genes that are expressed and
   have no multigene components, and the set of PacBio genes/transcripts that correpsond to our
   filtered gene set. Now I want to get the length of each feature and merge in.
   From here I will remove any hit that has <90% of the length of the feature, and then also
   flag if length of hit is 90%, 95% or 99% of original feature length */

data frag_length;
   set mm10.mm10_exon_fragment_flagged;
   feature_length=fragment_end-fragment_start;
   keep fragment_id feature_length;
   rename fragment_id=feature_id;
run;

data junc_length;
   set evspl.splicing_events_annot_refseq;
   feature_length=event_size;
   keep event_id feature_length;
   rename event_id=feature_id;
run;

data feat_length;
   set junc_length frag_length;
run;

proc sort data=feat_length;
   by feature_id;
proc sort data=blast_w_pb2refseq;
   by feature_id;
run;

data blast_w_feat_length;
   merge blast_w_pb2refseq (in=in1) feat_length (in=in2);
   by feature_id;
   if in1 and in2;
run;

/* Make permenant -- next we will filter on length, and flag hits if good/bad/etc. */

data event.pacbio_blast_hits_to_parse;
   set blast_w_feat_length;
run;







