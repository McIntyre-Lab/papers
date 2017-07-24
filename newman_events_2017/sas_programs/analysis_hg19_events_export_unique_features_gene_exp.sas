
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname splicing '/mnt/store/splice';
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';


/* Export counts for making plots to visualize the distribution of data post gene-level filtering
   I want to explore the following:
   Specifity*Feature type*Detection (post filtering!) -- thinking about a stacked barplot here!
   Features binned by #num xscripts
   e.g. 1 xs, 2 xs, 3 xs, 4 xs, 5xs, 6-10xs, 11-20xs, >20xs
*/

/* for each feature, get the type, detection (post-filtering), specificity and number of transcripts */

data feat_specificity;
   set event.hg19_flag_feat_specificity_exp;
   keep feature_id xscripts_per_feature flag_feature_unique flag_feature_common flag_feature_constitutive flag_multigene;
run;

data feat_detection;
  set event.hg19_features_w_annotations;
  where num_transcripts > 0;
  keep feature_id flag_cd4_on flag_cd8_on flag_cd19_on feature_type flag_singleton;
run;

data frag_to_keep;
   set event.hg19_flagged_fragment_length;
   where flag_fragment_lt_min_bp=0;
   keep fragment_id; rename fragment_id=feature_id;
run;

data junc_to_keep;
   set event.hg19_flagged_splicing_length;
   where flag_event_short=0;
   keep event_id; rename event_id=feature_id;
run;

data feat_to_keep;
   set junc_to_keep frag_to_keep;
run;



proc sort data=feat_specificity;
   by feature_id;
proc sort data=feat_detection;
   by feature_id;
proc sort data=feat_to_keep;
  by feature_id;
run;

data feat_data_for_plots short_feat;
  merge feat_to_keep (in=in1) feat_detection (in=in2) feat_specificity (in=in3);
  by feature_id;
  if in1 and in2 and in3 then output feat_data_for_plots;
  else if in2 and in3 then output short_feat;
run;

* Export;
proc export data=feat_data_for_plots outfile="!MCLAB/event_analysis/analysis_output/event_analysis_hg19_features_retained_with_info_gene_filtered_nomulti.csv" dbms=csv replace;
run;

data event.hg19_feature_data_for_plots_exp;
   set feat_data_for_plots;
run;

