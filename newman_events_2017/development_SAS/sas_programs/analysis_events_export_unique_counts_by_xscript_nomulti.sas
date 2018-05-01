ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count the number of unique features by type for genes and transcripts
   This is post length filtering. Need detection, but can use from above 

   Will have: gene/xscriptID, num_features, num_uniq_features, num_detected, num_uniq_features_detected
                              num_events
                              num_fragments*flag_singleton
                              num_fusions*flag_singleton
*/

data feat2xs2gene;
  set event.feature2xs2gene_nomulti;
run;

proc sort data=feat2xs2gene;
   by feature_id;
proc sort data=event.feature_data_for_plots;
   by feature_id;
run;

data feat_data_w_xs;
  merge feat2xs2gene (in=in1) event.feature_data_for_plots (in=in2);
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
  set event.feature2xs2gene_nomulti;
  keep transcript_id;
run;

data gene_list;
    set event.feature2xs2gene_nomulti;
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
   where flag_feature_on=1 and feature_type ne "fusion";
   var flag_feature_unique;
   output out=num_uniq_feat_dtct sum=num_unique_features_dtct;
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
   where flag_feature_on=1 and feature_type="splicing";
   var flag_feature_unique;
   output out=num_uniq_event_dtct sum=num_uniq_events_dtct;
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
   where flag_feature_on=1 and feature_type="fragment" and flag_singleton=1;
   var flag_feature_unique;
   output out=num_uniq_frag_sngl_dtct sum=num_uniq_frag_singletons_dtct;
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
   where flag_feature_on=1 and feature_type="fragment" and flag_singleton=0;
   var flag_feature_unique;
   output out=num_uniq_frag_nosngl_dtct sum=num_uniq_frag_nonsingl_dtct;
run;

* count number of unique fusion singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where feature_type="fusion" and flag_singleton=1;
   var flag_feature_unique;
   output out=num_uniq_fus_sngl sum=num_uniq_fus_singletons;
run;

* count number of detected unique fusion singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_feature_on=1 and feature_type="fusion" and flag_singleton=1;
   var flag_feature_unique;
   output out=num_uniq_fus_sngl_dtct sum=num_uniq_fus_singletons_dtct;
run;

* count number of unique fusion non-singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where feature_type="fusion" and flag_singleton=0;
   var flag_feature_unique;
   output out=num_uniq_fus_nosngl sum=num_uniq_fus_nonsingl;
run;

* count number of detected unique fusion non-singletons;
proc means data=&dataIn noprint;
   by &parentType._id;
   where flag_feature_on=1 and feature_type="fusion" and flag_singleton=0;
   var flag_feature_unique;
   output out=num_uniq_fus_nosngl_dtct sum=num_uniq_fus_nonsingl_dtct;
run;

* merge all together;
proc sort data=num_uniq_feat (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_feat_dtct (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_event (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_event_dtct (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_sngl (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_sngl_dtct (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_nosngl (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_frag_nosngl_dtct (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_fus_sngl (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_fus_sngl_dtct (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_fus_nosngl (drop=_TYPE_ _FREQ_);
  by &parentType._id;
proc sort data=num_uniq_fus_nosngl_dtct (drop=_TYPE_ _FREQ_);
  by &parentType._id;
run;

data &dataOut.;
  merge num_uniq_feat num_uniq_feat_dtct num_uniq_event num_uniq_event_dtct
        num_uniq_frag_sngl num_uniq_frag_sngl_dtct num_uniq_frag_nosngl
        num_uniq_frag_nosngl_dtct num_uniq_fus_sngl num_uniq_fus_sngl_dtct
        num_uniq_fus_nosngl num_uniq_fus_nosngl_dtct;
  by &parentType._id;
run;


proc sort data=&dataOut.;
   by &parentType._id;
proc sort data=&parentType._list;
   by &parentType._id;
run;

data &dataOut._2;
  merge &parentType._list (in=in1) &dataOut.;
  by &parentType._id;
  if in1;
run;

data event.&dataOut.;
   set &dataOut._2;
   array change _numeric_;
        do over change;
            if change=. then change=0;
        end;
run ;

%mend;

%count_features(feat_data_w_xs,transcript,uniq_feature_counts_by_xscript);
%count_features(feat2gene,gene,uniq_feature_counts_by_gene);

/* Export data */


proc export data=event.uniq_feature_counts_by_xscript
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_uniq_features_by_xscript_nomulti.csv" dbms=csv replace;
run;


proc export data=event.uniq_feature_counts_by_gene
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_uniq_features_by_gene_nomulti.csv" dbms=csv replace;
run;

