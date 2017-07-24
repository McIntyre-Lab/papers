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

* Fusions;

data fus2xscript;
   set mm10.mm10_si_fusions_unique_flagged;
   if num_exons_per_fusion=1 then flag_singleton=1;
   else if num_exons_per_fusion > 1 then flag_singleton=0;
   else delete;
   keep fusion_id flag_multigene flag_singleton num_xscripts_per_fusion;
run;

proc sort data=event.fusions_on_apn_gt0;
  by fusion_id;
proc sort data=fus2xscript;
  by fusion_id;
run;

data fus_w_xs;
   merge event.fusions_on_apn_gt0 (in=in1) fus2xscript (in=in2);
   by fusion_id;
   if in1 and in2;
run;

* Fragments ;

data frag2xscript;
    set mm10.mm10_refseq_exon_fragment_info;
    if num_genes_per_fragment > 1 then flag_multigene=1;
    else flag_multigene=0;
    keep fragment_id fusion_id num_xscripts_per_fragment flag_multigene;
run;

data fus_multi;
   set fus_w_xs;
   where flag_singleton=1;
   keep fusion_id;
run;

proc sort data=frag2xscript ;
   by fusion_id;
proc sort data=fus_multi;
   by fusion_id;
run;

data frag2xscript2;
   merge frag2xscript (in=in1) fus_multi (in=in2);
   by fusion_id;
   if in2 then flag_singleton=1;
   else flag_singleton=0;
   drop fusion_id;
run;

proc sort data=frag2xscript2;
   by fragment_id;
proc sort data=event.fragments_on_apn_gt0;
  by fragment_id;
run;

data frag_w_xs;
   merge event.fragments_on_apn_gt0 (in=in1) frag2xscript2 (in=in2);
   by fragment_id;
   if in1 and in2;
run;

* Splicing events;

data event2xscript;
   set evspl.splicing_events_annot_refseq;
   keep event_id num_transcripts;
run;

proc sort data=event2xscript;
  by event_id;
proc sort data=event.splicing_on_apn_gt0;
   by event_id;
run;

data event_w_xs;
  merge event.splicing_on_apn_gt0 (in=in1) event2xscript (in=in2);
  by event_id;
  if in1 and in2;
run;

/* Stack features and export so I can make some plots */

data fus_w_xs2;
   set fus_w_xs;
   length feature_type $8.;
   feature_type="fusion";
   flag_alt_acceptor=0;
   flag_alt_donor=0;
   flag_exonskip=0;
   flag_intron_retention=0;
   flag_junction_annotated=0;
   rename flag_fusion_on=flag_feature_on
          fusion_id=feature_id
          num_xscripts_per_fusion=num_transcripts;
run;

data frag_w_xs2;
   set frag_w_xs;
   length feature_type $8.;
   feature_type="fragment";
   flag_alt_acceptor=0;
   flag_alt_donor=0;
   flag_exonskip=0;
   flag_intron_retention=0;
   flag_junction_annotated=0;
   rename flag_fragment_on=flag_feature_on
          fragment_id=feature_id
          num_xscripts_per_fragment=num_transcripts;
run;

data event_w_xs2;
   set event_w_xs;
   length feature_type $8.;
   feature_type="splicing";
   flag_multigene=0;
   flag_singleton=0;
   rename flag_splicing_on=flag_feature_on
          event_id=feature_id;
run;

data event.features_w_annotations;
  set event_w_xs2 frag_w_xs2 fus_w_xs2;
run;

/* Export for plots */

proc export data=event.features_w_annotations 
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_features_w_annotations.csv" dbms=csv replace;
run;




