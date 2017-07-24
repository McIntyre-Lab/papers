ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname splicing '/mnt/store/splice';
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';


/* Export counts for making plots
   Add in the number of transcripts per feature
   For fusions and fragment, indicate if multigene and multiexon (in the case, if the REGION is multiexon)
*/

* Fragments ;

data frag2xscript;
    set hg19.hg19_aceview_exon_fragment_info;
    if num_genes_per_fragment > 1 then flag_multigene=1;
    else flag_multigene=0;
    keep fragment_id fusion_id num_xscripts_per_fragment flag_multigene;
run;

data fus_multi;
   set hg19.unique_info_fusions_si;
   if num_exons=1;
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

*update fragment IDs;

data frag_on;
   length chr $4.;
   set event.hg19_flag_fragment_on_apn5;
   chr=scan(fragment_id,1,":");
   fragment_start=scan(fragment_id,2,":") + 0;
   fragment_end=scan(fragment_id,3,":") + 0;
   rename fragment_id=fragment_coord;
run;

data frag_coord;
  set hg19.hg19_aceview_exon_Fragment_info;
  keep fragment_id chr fragment_start fragment_end;
run;

proc sort data=frag_on;
   by chr fragment_start fragment_end;
proc sort data=frag_coord;
   by chr fragment_start fragment_end;
run;

data frag_on2;
  merge frag_on (in=in1) frag_coord (in=in2);
  by chr fragment_start fragment_end;
  if in1 and in2;
run;


proc sort data=frag2xscript2;
   by fragment_id;
proc sort data=frag_on2;
  by fragment_id;
run;


data frag_w_xs;
   merge frag_on2 (in=in1) frag2xscript2 (in=in2);
   by fragment_id;
   if in1 and in2;
   drop chr fragment_start fragment_end fragment_coord;
run;

* Splicing events;

data event2xscript;
   set splicing.splicing_events_annotations;
   if num_transcripts > 0 then flag_junction_annotated=1;
   keep event_id num_transcripts flag_junction_annotated flag_intron_retention
        flag_exonskip flag_alt_donor flag_alt_acceptor;
run;

proc sort data=event2xscript;
  by event_id;
proc sort data=event.t1d_flag_splicing_on_apn5;
   by event_id;
run;

data event_w_xs;
  merge event.t1d_flag_splicing_on_apn5 (in=in1) event2xscript (in=in2);
  by event_id;
  if in1 and in2;
run;

/* Stack features and export so I can make some plots */

data frag_w_xs2;
   set frag_w_xs;
   length feature_type $8.;
   feature_type="fragment";
   flag_alt_acceptor=0;
   flag_alt_donor=0;
   flag_exonskip=0;
   flag_intron_retention=0;
   flag_junction_annotated=0;
   drop flag_fragment_all_on_ge5;
   rename fragment_id=feature_id
          num_xscripts_per_fragment=num_transcripts;
run;

data event_w_xs2;
   set event_w_xs;
   length feature_type $8.;
   feature_type="splicing";
   flag_multigene=0;
   flag_singleton=0;
   drop flag_splicing_all_ge5;
   rename event_id=feature_id;
run;

data event.hg19_features_w_annotations;
  set event_w_xs2 frag_w_xs2;
run;

/* Export for plots */

proc export data=event.hg19_features_w_annotations 
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_hg19_features_w_annotations.csv" dbms=csv replace;
run;




