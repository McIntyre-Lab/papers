ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count the number of unique features by type for genes and transcripts, post gene-level filtering
   This is post length filtering. Need detection, but can use from above 

   Will have: gene/xscriptID, num_features, num_uniq_features, num_detected, num_uniq_features_detected
                              num_events
                              num_fragments*flag_singleton
                              num_fusions*flag_singleton
*/

data feat2xs2gene;
  set event.hg19_feature2xs2gene_exp_only2;
run;

proc sort data=feat2xs2gene;
   by feature_id;
proc sort data=event.hg19_feature_data_for_plots_exp;
   by feature_id;
run;

data feat_data_w_xs;
  merge feat2xs2gene (in=in1) event.hg19_feature_data_for_plots_exp (in=in2);
  by feature_id;
  if in1 and in2;
run;

data feat2gene;
   set feat_data_w_xs;
   drop transcript_id;
run;

proc sort data=feat2gene nodup;
  by gene_id feature_id;
run;

data transcript_list;
  set event.hg19_feature2xs2gene_exp_only2;
  keep transcript_id;
run;

data gene_list;
    set event.hg19_feature2xs2gene_exp_only2;
  keep gene_id;
run;

proc sort data=transcript_list nodup;
   by transcript_id;
proc sort data=gene_list nodup;
   by gene_id;
run;

%macro count_features(dataIn,parentType,dataOut);

proc sort data=&dataIn.;
   by &parentType._id;
run;

* count number of unique features;
proc means data=&dataIn noprint;
   by &parentType._id;
   where feature_type ne "fusion";
   var flag_feature_unique;
   output out=num_uniq_feat sum=num_unique_features;
run;

* count number of detected unique features;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd4_on=1 and feature_type ne "fusion";
   var flag_feature_unique;
   output out=num_uniq_feat_dtct_cd4 sum=num_unique_features_dtct_cd4;
run;


* count number of detected unique features;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd8_on=1 and feature_type ne "fusion";
   var flag_feature_unique;
   output out=num_uniq_feat_dtct_cd8 sum=num_unique_features_dtct_cd8;
run;



* count number of detected unique features;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd19_on=1 and feature_type ne "fusion";
   var flag_feature_unique;
   output out=num_uniq_feat_dtct_cd19 sum=num_unique_features_dtct_cd19;
run;


* count number of unique events;
proc means data=&dataIn noprint;
   by &parentType._id;
   where feature_type="splicing";
   var flag_feature_unique;
   output out=num_uniq_event sum=num_uniq_events;
run;

* count number of detected unique events;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd4_on=1 and feature_type="splicing";
   var flag_feature_unique;
   output out=num_uniq_event_dtct_cd4 sum=num_uniq_events_dtct_cd4;
run;

* count number of detected unique events;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd8_on=1 and feature_type="splicing";
   var flag_feature_unique;
   output out=num_uniq_event_dtct_cd8 sum=num_uniq_events_dtct_cd8;
run;

* count number of detected unique events;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd19_on=1 and feature_type="splicing";
   var flag_feature_unique;
   output out=num_uniq_event_dtct_cd19 sum=num_uniq_events_dtct_cd19;
run;

* count number of unique fragment singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where feature_type="fragment" and flag_singleton=1;
   var flag_feature_unique;
   output out=num_uniq_frag_sngl sum=num_uniq_frag_singletons;
run;

* count number of detected unique fragment singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd4_on=1 and feature_type="fragment" and flag_singleton=1;
   var flag_feature_unique;
   output out=num_uniq_frag_sngl_dtct_cd4 sum=num_uniq_frag_sngltns_dtct_cd4;
run;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd8_on=1 and feature_type="fragment" and flag_singleton=1;
   var flag_feature_unique;
   output out=num_uniq_frag_sngl_dtct_cd8 sum=num_uniq_frag_sngltns_dtct_cd8;
run;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd19_on=1 and feature_type="fragment" and flag_singleton=1;
   var flag_feature_unique;
   output out=num_uniq_frag_sngl_dtct_cd19 sum=num_uniq_frag_sngltns_dtct_cd19;
run;

* count number of unique fragment non-singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where feature_type="fragment" and flag_singleton=0;
   var flag_feature_unique;
   output out=num_uniq_frag_nosngl sum=num_uniq_frag_nonsingl;
run;

* count number of detected unique fragment non-singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd4_on=1 and feature_type="fragment" and flag_singleton=0;
   var flag_feature_unique;
   output out=num_uniq_frag_nosngl_dtct_cd4 sum=num_uniq_frag_nonsingl_dtct_cd4;
run;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd8_on=1 and feature_type="fragment" and flag_singleton=0;
   var flag_feature_unique;
   output out=num_uniq_frag_nosngl_dtct_cd8 sum=num_uniq_frag_nonsingl_dtct_cd8;
run;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_cd19_on=1 and feature_type="fragment" and flag_singleton=0;
   var flag_feature_unique;
   output out=num_uniq_frag_nosngl_dtct_cd19 sum=num_uniq_frag_nonsingl_dtct_cd19;
run;

* merge all together;
proc sort data=num_uniq_feat (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_feat_dtct_cd4 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_feat_dtct_cd8 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_feat_dtct_cd19 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_event (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_event_dtct_cd4 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_event_dtct_cd8 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_event_dtct_cd19 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_sngl (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_sngl_dtct_cd4 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_sngl_dtct_cd8 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_sngl_dtct_cd19 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_nosngl (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_nosngl_dtct_cd4 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_nosngl_dtct_cd8 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_nosngl_dtct_cd19 (drop=_TYPE_ _FREQ_);
  by &parentType._id;
run;

data merged_data;
  merge num_uniq_feat num_uniq_feat_dtct_cd4 num_uniq_feat_dtct_cd8 num_uniq_feat_dtct_cd19
        num_uniq_event num_uniq_event_dtct_cd4 num_uniq_event_dtct_cd8 num_uniq_event_dtct_cd19
        num_uniq_frag_sngl num_uniq_frag_sngl_dtct_cd4 num_uniq_frag_sngl_dtct_cd8 num_uniq_frag_sngl_dtct_cd19
        num_uniq_frag_nosngl num_uniq_frag_nosngl_dtct_cd4 num_uniq_frag_nosngl_dtct_cd8
         num_uniq_frag_nosngl_dtct_cd19;
  by &parentType._id;
run;


proc sort data=merged_data;
   by &parentType._id;
proc sort data=&parentType._list;
   by &parentType._id;
run;

data merged_data_2;
  merge &parentType._list (in=in1) merged_data;
  by &parentType._id;
  if in1;
run;


data &dataOut.;
   set merged_data_2;
   array change _numeric_;
        do over change;
            if change=. then change=0;
        end;
run ;

%mend;


%count_features(feat_data_w_xs,transcript,hg19_feat_cnt_by_xs_gene_exp);
%count_features(feat2gene,gene,hg19_feat_cnt_by_gene_exp);

/* Make permenant */

data event.hg19_feat_cnt_by_xs_gene_exp;
   set hg19_feat_cnt_by_xs_gene_exp;
run;

data event.hg19_feat_cnt_by_gene_exp;
   set hg19_feat_cnt_by_gene_exp;
run;




/* Export data */

proc export data=event.hg19_feat_cnt_by_xs_gene_exp
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_hg19_uniq_features_by_xscript_gene_filtered_nomulti.csv" dbms=csv replace;
run;

proc export data=event.hg19_feat_cnt_by_gene_exp
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_hg19_uniq_features_by_gene_gene_filtered_nomulti.csv" dbms=csv replace;
run;

