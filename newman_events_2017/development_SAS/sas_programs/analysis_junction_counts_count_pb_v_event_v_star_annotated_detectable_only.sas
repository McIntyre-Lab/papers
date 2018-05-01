/* Counts needed:
(5) Annotated junctions:
	Restrict to ONLY PB junctions, and compare to PB junctions
		Also to genes with PB transcripts
	Count junctions detected by each method
	STAR/events vs PB
	STAR vs Events
*/
ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Using PacBio as a reference, compare the number of annotated junctions identified with EA or STAR
   against the annotated PacBio junctions

   (1) Annotated EA vs PacBio
   (2) Annotated STAR vs PacBio
   (3) Annotated EA vs STAR */

/* Create an index number for junctions */

data junc_index;
  set event.catalog_pacbio_star_junctions;
  length junc_index $12.;
  junc_index=catt('junc',_n_);
  if flag_apn_NSC1_gt0=1 or flag_apn_NSC2_gt0=1 then flag_NPC_any_apn0=1; else flag_NPC_any_apn0=0;
  if flag_STAR_depth_NSC1_gt0=1 or flag_STAR_depth_NSC2_gt0=1 then flag_NPC_any_depth0=1; else  flag_NPC_any_depth0=0;
run;

data pb_annot_junc;
   set junc_index;
   if flag_in_pacbio=1 and flag_junction_annotated = 1;
   keep junc_index;
run;

data cat_annot;
   set junc_index;
   where flag_NPC_any_apn0=1;
   if flag_in_catalog=1 and flag_junction_annotated=1;
   keep junc_index;
run;

data star_annot;
   set junc_index;
  where flag_NPC_any_depth0=1;
   if flag_in_star=1  and flag_junction_annotated=1;
   keep junc_index;
run;

data detection;
   set junc_index;
   if flag_apn_NSC1_gt0=1 and flag_apn_NSC2_gt0=1 then flag_NPC_both_apn0=1; else flag_NPC_both_apn0=0;
   if flag_apn_NSC1_ge5=1 and flag_apn_NSC2_ge5=1 then flag_NPC_both_apn5=1; else flag_NPC_both_apn5=0;
   if flag_STAR_depth_NSC1_gt0=1 and flag_STAR_depth_NSC2_gt0=1 then flag_NPC_both_depth0=1; else flag_NPC_both_depth0=0;
   if flag_STAR_depth_NSC1_ge5=1 and flag_STAR_depth_NSC2_ge5=1 then flag_NPC_both_depth5=1; else flag_NPC_both_depth5=0;
    keep junc_index flag_NPC_both_apn0 flag_NPC_both_apn5
                    flag_NPC_both_depth0 flag_NPC_both_depth5
      chr donor_stop acceptor_start strand event_id;
run;

data juncs_pb_regions;
   set event.catalog_pacbio_star_junc_pb;
   keep chr donor_stop acceptor_start strand event_id flag_in_pb_coord;
run;

proc sort data=detection;
   by chr strand donor_stop acceptor_start event_id;
proc sort data=juncs_pb_regions;
   by chr strand donor_stop acceptor_start event_id;
run;

data detection2;
  merge detection (in=in1) juncs_pb_regions (in=in2);
  by chr strand donor_stop acceptor_start event_id;
  if in2 then flag_gene_has_pb=1; else flag_gene_has_pb=0;
  if in1 then output;
run;

/* Compare STAR to PacBio -- annotated junctions */
/***** Annotated : STAR vs Pacbio *****/

%include "!MCLAB/event_analysis/sas_programs/analysis_junc_counts_star_vs_pacbio_annotated_detectable_only.sas";

/***** Annotated : EA vs Pacbio *****/
%include "!MCLAB/event_analysis/sas_programs/analysis_junc_counts_events_vs_pacbio_annotated_detectable_only.sas";
/***** EA vs STAR: Annot EA vs Annot STAR
Look at: all possible, PB junctions only, Junction in genes with PB transcripts

*****/
%include "!MCLAB/event_analysis/sas_programs/analysis_junc_counts_events_vs_star_annotated_detectable_only.sas";








