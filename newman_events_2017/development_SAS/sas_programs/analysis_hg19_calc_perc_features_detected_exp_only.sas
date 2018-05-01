ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Creating the list of transcripts (from expressed genes) that are possibly expressed, based on the exon/fragment/splicing expressed data. A threshold of the proportion of unique features that are detected is set (user-defined). For the analysis here, I am using 4 values: 25%, 50%, 75% and 100%. For transcripts with no unique features, these automatically are put into the list of possibly expressed transcripts (cannot resolve)

Flag any transcript with at least 1 unique feature, and at least 1 unique feature detected, and calculate the proportion of unique features detected of all assigned unique features.
*/



(1) feature specificity (to calculate #uniq and %uniq dtct)
	event.hg19_flag_feat_specificity_exp;
(2) feature-to-xs index
	event.feature2xs2gene_exp_only2
(3) detection flags
	event.hg19_features_w_annotations




data feat2xs2gene;
  set event.hg19_feature2xs2gene_exp_only2;
run;

data uniq_features;
   set event.hg19_flag_feat_specificity_exp;
   where flag_feature_unique=1;
   keep feature_id;
run;

data flag_on;
   set event.hg19_features_w_annotations;
   keep feature_id flag_cd4_on flag_cd8_on flag_cd19_on;
run;


proc sort data=feat2xs2gene;
   by feature_id;
proc sort data=uniq_features;
   by feature_id;
proc sort data=flag_on;
   by feature_id;
run;

data feat_data_w_xs;
  merge feat2xs2gene (in=in1) uniq_features (in=in2) flag_on (in=in3);
  by feature_id;
  if in2 then flag_feature_unique=1; else flag_feature_unique=0;
  if not in3 then do;
     flag_cd4_on=0; flag_cd8_on=0; flag_cd19_on=0; end;
  flag_feature=1; *set this here so I can just sum this below;
  if in1 then output;
run;

proc sort data=feat_data_w_xs;
   by transcript_id;
run;

* count number of features;
proc means data=feat_data_w_xs noprint;
   by transcript_id;
   var flag_feature;
   output out=num_feat sum=num_total_features;
run;

* count number of detected features;
proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_cd4_on=1;
   var flag_feature;
   output out=num_feat_dtct_cd4 sum=num_features_dtct_cd4;
run;

proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_cd8_on=1;
   var flag_feature;
   output out=num_feat_dtct_cd8 sum=num_features_dtct_cd8;
run;

proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_cd19_on=1;
   var flag_feature;
   output out=num_feat_dtct_cd19 sum=num_features_dtct_cd19;
run;

* count number of unique features;
proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_feature_unique=1;
   var flag_feature;
   output out=num_uniq_feat sum=num_total_uniq_features;
run;

* count number of detected unique features;
proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_cd4_on=1 and  flag_feature_unique=1;
   var flag_feature;
   output out=num_uniq_feat_dtct_cd4 sum=num_uniq_features_dtct_cd4;
run;

proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_cd8_on=1 and  flag_feature_unique=1;
   var flag_feature;
   output out=num_uniq_feat_dtct_cd8 sum=num_uniq_features_dtct_cd8;
run;

proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_cd19_on=1 and  flag_feature_unique=1;
   var flag_feature;
   output out=num_uniq_feat_dtct_cd19 sum=num_uniq_features_dtct_cd19;
run;


* merge all together;
proc sort data=num_feat (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_feat_dtct_cd4 (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_feat_dtct_Cd8 (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_feat_dtct_cd19 (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_uniq_feat (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_uniq_feat_dtct_cd4 (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_uniq_feat_dtct_Cd8 (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_uniq_feat_dtct_cd19 (drop=_TYPE_ _FREQ_);
  by transcript_id;
run;

data num_features_by_xs;
  merge num_feat num_feat_dtct_cd4 num_feat_dtct_cd8  num_feat_dtct_cd19
      num_uniq_feat num_uniq_feat_dtct_cd4 num_uniq_feat_dtct_cd8  num_uniq_feat_dtct_cd19;
  by transcript_id;
run;

data num_features_by_xs3;
   set num_features_by_xs;
   array change _numeric_;
        do over change;
            if change=. then change=0;
        end;
run ;

/* Calculate the proportion of features detected, and the proportion of unique features detected */

data calc_prop_dtct;
  set num_features_by_xs3;
  perc_features_dtct_cd4=num_features_dtct_cd4/num_total_features;
  perc_features_dtct_cd8=num_features_dtct_cd8/num_total_features;
  perc_features_dtct_cd19=num_features_dtct_cd19/num_total_features;
  if num_total_uniq_features ne 0 then do;
        perc_unique_features_dtct_cd4=num_uniq_features_dtct_cd4 / num_total_uniq_features;
        perc_unique_features_dtct_cd8=num_uniq_features_dtct_cd8 / num_total_uniq_features;
        perc_unique_features_dtct_cd19=num_uniq_features_dtct_cd19 / num_total_uniq_features;
  end;
run;

* Flag transcripts if they have uniq, unique on, etc.;

data flag_xscripts;
  set calc_prop_dtct;
  if num_total_features > 0 then flag_xscript_has_features=1;
     else flag_xscript_has_features=0;
  if sum(num_features_dtct_cd4,num_features_dtct_cd8,num_features_dtct_cd19) > 0 then flag_xscript_has_dtct_features=1;
     else flag_xscript_has_dtct_features=0;
  if num_total_uniq_features > 0 then flag_xscript_has_unique=1;
     else flag_xscript_has_unique=0;
  if sum(num_uniq_features_dtct_cd4,num_uniq_features_dtct_cd8,num_uniq_features_dtct_cd19) > 0 then flag_xscript_has_unique_dtct=1;
     else flag_xscript_has_unique_dtct=0;
run;

/* Make permenant */

data event.hg19_flag_xscripts_w_unique;
  set flag_xscripts;
run;


/* Checks: number of xscripts with features, features detected, unique features */


/*
proc freq data=flag_xscripts;
   tables flag_xscript_has_features flag_xscript_has_dtct_features
          flag_xscript_has_unique flag_xscript_has_unique_dtct;
run;
*/

/*
 
   flag_xscript_                             Cumulative    Cumulative
    has_features    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               1       73535      100.00         73535       100.00


   flag_xscript_
       has_dtct_                             Cumulative    Cumulative
        features    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               1       73535      100.00         73535       100.00


   flag_xscript_                             Cumulative    Cumulative
      has_unique    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       34900       47.46         34900        47.46
               1       38635       52.54         73535       100.00


   flag_xscript_                             Cumulative    Cumulative
 has_unique_dtct    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       51216       69.65         51216        69.65
               1       22319       30.35         73535       100.00

*/



