ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Exclude events (junctions, exon fragments) and transcripts from genes with any multigene components 

Create new set of annotations, excluding (1) genes with multigene components (2) ambiguous junctions
*/


ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Export counts for making plots
   Add in the number of transcripts per feature
   For fusions and fragment, indicate if multigene and multiexon (in the case, if the REGION is multiexon)
*/


data genes_to_keep;
  set event.genes_flag_multigene;
  where flag_multigene_region=0;
  keep gene_id;
run;

data ambig_junc;
   set evspl.flag_event_redundant_seq;
   where flag_redundant_sequence=1;
   keep event_id;
   rename event_id=feature_id;
run;

data feat2xs2gene;
   set event.feature2xs2gene;
run;

proc sort data=feat2xs2gene;
   by feature_id;
proc sort data=ambig_junc;
   by feature_id;
run;

data feat2xs2gene_2;
  merge feat2xs2gene (in=in1) ambig_junc (in=in2);
  by feature_id;
  if in2 then delete;
run;

proc sort data=feat2xs2gene_2;
   by gene_id;
proc sort data=genes_to_keep;
  by gene_id;
run;

data feat2xs2gene_3 dropped;
  merge feat2xs2gene_2 (in=in1) genes_to_keep (in=in2);
  by gene_id;
  if in1 and in2 then output feat2xs2gene_3;
  else if in1 then output dropped;
run;

data feat2keep;
  set feat2xs2gene_3;
  keep feature_id;
run;

data unannot_feat;
   set evspl.splicing_events_annot_refseq;
   keep gene_id event_id;
   rename event_id=feature_id;
run;

proc sort data=unannot_feat;
  by feature_id;
run;

data unannot_feat2;
  merge unannot_feat (in=in1) ambig_junc (in=in2);
  by feature_id;
  if in2 then delete;
run;

proc sort data=unannot_feat2;
  by gene_id;
proc sort data=genes_to_keep;
   by gene_id;
run;

data unannot_feat2keep;
   merge unannot_feat2 (in=in1) genes_to_keep (in=in2);
   by gene_id;
   if in1 and in2;
run;

data unannot_feat2keep2;
   set unannot_feat2keep;
   keep feature_id;
run;

data feat2keep2;
   set feat2keep unannot_feat2keep2;
run;

data feat_annot;
   set event.features_w_annotations;
run;

proc sort data=feat_annot;
  by feature_id;
proc sort data=feat2keep2 nodup;
  by feature_id;
run;

data feat_annot2 dropped;
  merge feat_annot (in=in1) feat2keep2 (in=in2);
  by feature_id;
  if in1 and in2 then output feat_annot2;
  else output dropped;
run;


/* Make permenant */

data event.feature2xs2gene_nomulti;
   set feat2xs2gene_3;
run;

data event.features_w_annotations_nomulti;
   set feat_annot2;
run;



/* Export for plots */

proc export data=event.features_w_annotations_nomulti 
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_features_w_annotations_nomulti.csv" dbms=csv replace;
run;


