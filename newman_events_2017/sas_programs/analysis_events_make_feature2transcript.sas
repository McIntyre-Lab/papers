ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* From the list of filtered features, create a list of features by transcript and gene */

* List of features after length filtering;
data feature_gt_min;
  set event.flagged_feature_short;
  where flag_feature_short=0;
  keep feature_id;
run;

* Fusion to transcript;
data fus2xs;
  set mm10.mm10_fusion_si_info_unique;
  length transcript_id2 $20.;
  if transcript_id = '' then delete;
  else do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      transcript_id2=scan(transcript_id,i,"|");
      keep fusion_id transcript_id2;
      output;
      end;
  rename fusion_id=feature_id transcript_id2=transcript_id;
run;


* Fragment to transcript;
data frag2xs;
  set mm10.mm10_exon_fragment_flagged;
  length transcript_id2 $20.;
  if transcript_id = '' then delete;
  else do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      transcript_id2=scan(transcript_id,i,"|");
      keep fragment_id transcript_id2;
      output;
      end;
  rename fragment_id=feature_id transcript_id2=transcript_id;
run;


* Splicing to transcript;
data event2xs;
  set evspl.splicing_events_annot_refseq;
  length transcript_id2 $20.;
  if transcript_id = '' then delete;
  else do i=1 by 1 while(scan(transcript_id,i,"|") ^= "");
      transcript_id2=scan(transcript_id,i,"|");
      keep event_id transcript_id2;
      output;
      end;
  rename event_id=feature_id transcript_id2=transcript_id;
run;


* Transcript to gene;
data xs2gene;
  set refseq.ref_transcript2geneid;
run;

* Stack features;
data feature2xs;
  length feature_id $450.;
  format feature_id $450.;
  set event2xs frag2xs fus2xs;
run;

proc sort data=xs2gene;
   by transcript_id;
proc sort data=feature2xs;
   by transcript_id;
run;

data feature2xs2gene;
   merge feature2xs (in=in1) xs2gene (in=in2);
   by transcript_id;
   if in1 and in2;
run;

/* Make permenant */

data event.feature2xs2gene;
  set feature2xs2gene;
run;

