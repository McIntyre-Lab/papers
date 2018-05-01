/* Summarizing some event analysis stuff:

(3) Calculate the distribution of transcripts by number of detected unique pieces
	(a) #uniq * transcript (features, fragments, events, fusions)
	(b) %uniq * transcript (features, fragments, events, fusions)

  For each transcript, calculate the number unique features detected by type, and the proportion of unique features of detected features:
    num_uniq_fragments
    num_uniq_fusions
    num_uniq_junctions
    num_uniq_frags_juncs (features: fragments + junctions)
    num_uniq_fus_juncs (features: fusions + junctions)
    num_dtct_fragments
    num_dtct_fusions
    num_dtct_junctions
    num_dtct_frags_juncs  (features: fragments + junctions)
    num_dtct_fus_juncs  (features: fusions + junctions)

  I am going to output this as a CSV and plot in Python
*/

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Get features to transcripts */

data xs_w_uniq;
  set event.flag_xscript_uniq_features;
  keep transcript_id;
run;

data frag2xs;
   set mm10.mm10_exon_fragment_flagged;
   length transcript_id2 $20.;
   length feature_type $8.;
   feature_type="fragment";
   if fragment_end-fragment_start < 27 then delete;
   do i=1 by 1 while(scan(transcript_id,i,"|") ^= " ");
       transcript_id2=scan(transcript_id,i,"|");
       keep fragment_id transcript_id2 feature_type;
       rename fragment_id=feature_id transcript_id2=transcript_id ;
       output;
       end;
run;

data fus_event2xs;
   set refseq.ref_fusions_junctions_stack;
   drop gene_id;
run;

data feature2xs;
   set fus_event2xs frag2xs;
run;

proc sort data=feature2xs;
   by transcript_id;
proc sort data=xs_w_uniq;
  by transcript_id;
run;

data xs_w_uniq_all_feats;
  merge xs_w_uniq (in=in1) feature2xs (in=in2);
  by transcript_id;
  if in1 and in2;
run;

/* Get features detected by commonality */

data detected_features;
  set event.dtct_features_by_class;
  if flag_feature_unique=1 then flag_feature_constitutive=0;
  keep feature_id flag_feature_unique flag_feature_common flag_feature_constitutive; 
run;

proc sort data=detected_features nodup;
  by feature_id;
proc sort data=xs_w_uniq_all_feats;
  by feature_id;
run;

data xs_w_uniq_dtct_feats;
   merge xs_w_uniq_all_feats (in=in1) detected_features (in=in2);
   by feature_id;
   if in1 and in2;
   flag_feature_detected=1;
run;

/* Count number of fragments by type */

proc sort data=xs_w_uniq_dtct_feats;
   by transcript_id;
proc means data=xs_w_uniq_dtct_feats noprint;
   by transcript_id;
   where feature_type="fragment";
   var flag_feature_unique flag_feature_common flag_feature_constitutive flag_feature_detected;
   output out=num_frag_by_type
          sum(flag_feature_unique)=num_unique_fragment
          sum(flag_feature_common)=num_common_fragment
          sum(flag_feature_constitutive)=num_constit_fragment
          sum(flag_feature_detected)=num_dtct_fragment;
run;

/* Count number of fusions by type */

proc means data=xs_w_uniq_dtct_feats noprint;
   by transcript_id;
   where feature_type="fusion";
   var flag_feature_unique flag_feature_common flag_feature_constitutive flag_feature_detected;
   output out=num_fus_by_type
          sum(flag_feature_unique)=num_unique_fusion
          sum(flag_feature_common)=num_common_fusion
          sum(flag_feature_constitutive)=num_constit_fusion
          sum(flag_feature_detected)=num_dtct_fusion;
run;

/* Count number of splicing by type */

proc means data=xs_w_uniq_dtct_feats noprint;
   by transcript_id;
   where feature_type="junction";
   var flag_feature_unique flag_feature_common flag_feature_constitutive flag_feature_detected;
   output out=num_junc_by_type
          sum(flag_feature_unique)=num_unique_junction
          sum(flag_feature_common)=num_common_junction
          sum(flag_feature_constitutive)=num_constit_junction
          sum(flag_feature_detected)=num_dtct_junction;
run;

proc sort data=num_frag_by_type;
  by transcript_id;
proc sort data=num_fus_by_type;
  by transcript_id;
proc sort data=num_junc_by_type;
  by transcript_id;
run;

data event.num_dtct_features_by_isoform;
   merge num_frag_by_type (in=in1) num_fus_by_type (in=in2) num_junc_by_type (in=in3);
   length isoform_group_feat1 $20.;
   length isoform_group_feat2 $20.;
   by transcript_id;
   if not in1 then do;
      num_unique_fragment=0;
      num_common_fragment=0;
      num_constit_fragment=0;
      num_dtct_fragment=0; end;
   if not in2 then do;
      num_unique_fusion=0;
      num_common_fusion=0;
      num_constit_fusion=0;
      num_dtct_fusion=0; end;
   if not in3 then do;
      num_unique_junction=0;
      num_common_junction=0;
      num_constit_junction=0;
      num_dtct_junction=0; end;
   num_unique_feature1=num_unique_fragment+num_unique_junction;
   num_common_feature1=num_common_fragment+num_common_junction;
   num_constit_feature1=num_constit_fragment+num_constit_junction;
   num_dtct_feature1=num_dtct_fragment+num_dtct_junction;
   if num_unique_feature1 > 0 then isoform_group_feat1="unique";
   else if num_common_feature1 > 0  then isoform_group_feat1="common";
   else isoform_group_feat1="constitutive";

   num_unique_feature2=num_unique_fusion+num_unique_junction;
   num_common_feature2=num_common_fusion+num_common_junction;
   num_constit_feature2=num_constit_fusion+num_constit_junction;
   num_dtct_feature2=num_dtct_fusion+num_dtct_junction;
   if num_unique_feature2 > 0 then isoform_group_feat2="unique";
   else if num_common_feature2 > 0  then isoform_group_feat2="common";
   else isoform_group_feat2="constitutive";
   drop _TYPE_ _FREQ_;
run;

* Export for python;

proc export data=event.num_dtct_features_by_isoform
    outfile="!MCLAB/event_analysis/analysis_output/detected_features_by_isoform.csv"
    dbms=csv replace;
run;





  
