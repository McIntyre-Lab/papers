ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname splicing '/mnt/store/splice';
libname eventloc '/mnt/store/event_sandbox/sas_data';

/* Flag hg19 junctions and fragments that are too short.
   For fragments this is a minimum length of 10bp
   For junctions this is the flag_event_short flag */

data frag2keep;
   set event.hg19_feature2xs2gene_exp_only;
   keep feature_id;
   rename feature_id=fragment_id;
run;

data gene2keep;
  set event.hg19_flag_gene_expressed;
  if flag_multigene=1 then delete;
  if flag_cd4_gene_on=1 or flag_cd8_gene_on=1 or flag_cd19_gene_on=1 then output;
  keep gene_id;
run; *27993 genes kept;


data frag_length;
  set hg19.hg19_aceview_exon_fragment_info;
  fragment_length=fragment_end-fragment_start;
  if fragment_length < 10 then flag_fragment_lt_min_bp=1;
   else flag_fragment_lt_min_bp=0;
  keep fragment_id flag_fragment_lt_min_bp;
run;

data junc_length;
  set splicing.splicing_events_annotations;
  if event_size < 10 then flag_splicing_lt_min_bp=1;
   else flag_splicing_lt_min_bp=0;
  keep event_id flag_splicing_lt_min_bp flag_event_short gene_id;
run;

 
proc sort data=gene2keep nodup;
  by gene_id;
proc sort data=junc_length;
  by gene_id;
run;

data junc_length2;
  merge junc_length (in=in1) gene2keep (in=in2);
  by gene_id;
  if in1 and in2;
  drop gene_id;
run;

 
proc sort data=frag2keep nodup;
  by fragment_id;
proc sort data=frag_length;
  by fragment_id;
run;

data frag_length2;
  merge frag2keep  (in=in1) frag_length (in=in2);
  by fragment_id;
  if in1 and in2;
run;

/* Remove events that are short from event-to-xscript-to-gene index */

data event2xs;
  set event.hg19_feature2xs2gene_exp_only;
run;

data frag_short;
  set frag_length2;
  where flag_fragment_lt_min_bp=1;
  keep fragment_id;
  rename fragment_id=feature_id;
run;


data junc_short;
  set junc_length2;
  where flag_event_short=1;
  keep event_id;
  rename event_id=feature_id;
run;

data feat_short;
  set junc_short frag_short;
run;

proc sort data=feat_short;
  by feature_id;
proc sort data=event2xs;
  by feature_id;
run;

data event2xs_filt;
  merge event2xs (in=in1) feat_short (in=in2);
  by feature_id;
  if in2 then delete;
run;


/* Make permenant */

data event.hg19_flagged_fragment_length;
   set frag_length2;
run;

data eventloc.hg19_flagged_splicing_length;
   set junc_length2;
run;

data event.hg19_feature2xs2gene_exp_only2;
   set event2xs_filt;
run;

