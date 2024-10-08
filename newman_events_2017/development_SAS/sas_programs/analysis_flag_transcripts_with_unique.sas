ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Creating the list of transcripts (from expressed genes) that are possibly expressed, based on the exon/fragment/splicing expressed data. A threshold of the proportion of unique features that are detected is set (user-defined). For the analysis here, I am using 4 values: 25%, 50%, 75% and 100%. For transcripts with no unique features, these automatically are put into the list of possibly expressed transcripts (cannot resolve)

Flag any transcript with at least 1 unique feature, and at least 1 unique feature detected, and calculate the proportion of unique features detected of all assigned unique features.
*/

/* In the last program I calculated the number of unique and detected unique per transcript. I want to also calculate the number of total features detected, in case there are retained transcripts with nothing detected */

data feat2xs2gene;
  set event.feature2xs2gene_exp_only;
run;

proc sort data=feat2xs2gene;
   by feature_id;
proc sort data=event.feature_data_for_plots_gene_exp;
   by feature_id;
run;

data feat_data_w_xs;
  merge feat2xs2gene (in=in1) event.feature_data_for_plots_gene_exp (in=in2);
  by feature_id;
  if in1 and in2;
  flag_feature=1; *set this here so I can just sum this below;
run;

proc sort data=feat_data_w_xs;
   by transcript_id;
run;

* count number of features;
proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where feature_type ne "fusion";
   var flag_feature;
   output out=num_feat sum=num_total_features;
run;

* count number of detected features;
proc means data=feat_data_w_xs noprint;
   by transcript_id;
   where flag_feature_on=1 and feature_type ne "fusion";
   var flag_feature;
   output out=num_feat_dtct sum=num_features_dtct;
run;

* merge all together;
proc sort data=num_feat (drop=_TYPE_ _FREQ_);
  by transcript_id;
proc sort data=num_feat_dtct (drop=_TYPE_ _FREQ_);
  by transcript_id;
run;

data num_features_by_xs;
  merge num_feat num_feat_dtct;
  by transcript_id;
run;

/* Merge in with unique feature counts */

data uniq_counts;
   set event.feature_cnt_by_xscript_gene_exp;
   keep transcript_id num_unique_features num_unique_features_dtct;
run;

proc sort data=num_features_by_xs;
  by transcript_id;
proc sort data=uniq_counts;
  by transcript_id;
run;

data num_features_by_xs2;
  merge num_features_by_xs uniq_counts (in=in2);
  by transcript_id;
  if in2;
run;

data num_features_by_xs3;
   set num_features_by_xs2;
   array change _numeric_;
        do over change;
            if change=. then change=0;
        end;
run ;

/* Calculate the proportion of features detected, and the proportion of unique features detected */

data calc_prop_dtct;
  set num_features_by_xs3;
  perc_features_dtct=num_features_dtct/num_total_features;
  if num_unique_features ne 0 then perc_unique_features_dtct=num_unique_features_dtct / num_unique_features;
run;

* Flag transcripts if they have uniq, unique on, etc.;

data flag_xscripts;
  set calc_prop_dtct;
  if num_total_features > 0 then flag_xscript_has_features=1;
     else flag_xscript_has_features=0;
  if num_features_dtct > 0 then flag_xscript_has_dtct_features=1;
     else flag_xscript_has_dtct_features=0;
  if num_unique_features > 0 then flag_xscript_has_unique=1;
     else flag_xscript_has_unique=0;
  if num_unique_features_dtct > 0 then flag_xscript_has_unique_dtct=1;
     else flag_xscript_has_unique_dtct=0;
run;

/* Make permenant */

data event.flag_xscripts_w_unique;
  set flag_xscripts;
run;



