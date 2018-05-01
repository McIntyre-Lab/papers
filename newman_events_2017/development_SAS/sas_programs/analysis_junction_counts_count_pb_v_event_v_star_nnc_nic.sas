/* Counts needed:
(1) PacBio NNC and NIC junctions
(2) Novel ID rate:
	Unannot EA vs STAR novel
	Unannot+Border EA vs STAR novel
(3) Redo counts, include ALL PB novel in PB-detected set
(4) Novel counts: NIC: STAR/EA vs PB
			STAR/EA: annotated as "unannot" in event catalog, vs NIC PB
                  NNC: STAR vs PB
			STAR: not in event catalog, vs NIC PB
(5) Annotated junctions:
	Restrict to ONLY PB junctions, and compare to PB junctions
		Also to genes with PB transcripts
	Count junctions detected by each method
	STAR/events vs PB
	STAR vs Events
*/

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/splicing/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Using PacBio as a reference, compare the number of novel/unannotated junctions identified with EA or STAR
   against the NNC/NIC PacBio junctions

   (1) Unannotated EA vs NIC PacBio
   (2) Unannotated STAR vs NIC PacBio
   (3) Novel STAR vs NNC PacBio */

/* Create an index number for junctions */

data junc_index;
  set event.catalog_pacbio_star_junctions;
  length junc_index $12.;
  junc_index=catt('junc',_n_);
run;

data pb_nic_junc pb_nnc_junc;
   set junc_index;
   if flag_in_pacbio=1 and flag_junction_annotated ne 1 then do;
   if flag_junction_annotated=0 then output pb_nic_junc;
   if flag_junction_annotated=. then output pb_nnc_junc;
   end;
   keep junc_index;
run;

data cat_unannot;
   set junc_index;
   if flag_in_catalog=1 and flag_junction_annotated=0;
   keep junc_index;
run;

data star_unannot star_novel;
   set junc_index;
   if flag_in_star=1 then do;
      if flag_junction_annotated=0 then output star_unannot;
      if flag_junction_annotated=. then output star_novel;
   end;
   keep junc_index;
run;

data detection;
   set junc_index;
   if flag_apn_NSC1_gt0=1 or flag_apn_NSC2_gt0=1 then flag_NPC_any_apn0=1; else flag_NPC_any_apn0=0;
   if flag_apn_NSC1_gt0=1 and flag_apn_NSC2_gt0=1 then flag_NPC_both_apn0=1; else flag_NPC_both_apn0=0;
   if flag_apn_NSC1_ge5=1 or flag_apn_NSC2_ge5=1 then flag_NPC_any_apn5=1; else flag_NPC_any_apn5=0;
   if flag_apn_NSC1_ge5=1 and flag_apn_NSC2_ge5=1 then flag_NPC_both_apn5=1; else flag_NPC_both_apn5=0;

   if flag_STAR_depth_NSC1_gt0=1 or flag_STAR_depth_NSC2_gt0=1 then flag_NPC_any_depth0=1; else flag_NPC_any_depth0=0;
   if flag_STAR_depth_NSC1_gt0=1 and flag_STAR_depth_NSC2_gt0=1 then flag_NPC_both_depth0=1; else flag_NPC_both_depth0=0;
   if flag_STAR_depth_NSC1_ge5=1 or flag_STAR_depth_NSC2_ge5=1 then flag_NPC_any_depth5=1; else flag_NPC_any_depth5=0;
   if flag_STAR_depth_NSC1_ge5=1 and flag_STAR_depth_NSC2_ge5=1 then flag_NPC_both_depth5=1; else flag_NPC_both_depth5=0;
    keep junc_index flag_NPC_any_apn0 flag_NPC_both_apn0 flag_NPC_any_apn5 flag_NPC_both_apn5
                    flag_NPC_any_depth0 flag_NPC_both_depth0 flag_NPC_any_depth5 flag_NPC_both_depth5
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


/***** NIC : STAR vs Pacbio *****/
%include "!MCLAB/event_analysis/sas_programs/analysis_junc_counts_pb_vs_star_nic.sas";


/***** NIC : EA vs Pacbio *****/
%include "!MCLAB/event_analysis/sas_programs/analysis_junc_counts_pb_vs_events_nic.sas";

/***** NNC : STAR vs Pacbio *****/
%include "!MCLAB/event_analysis/sas_programs/analysis_junc_counts_pb_vs_star_nnc.sas";

/***** EA vs STAR:
(1) Unannot EA vs Unannot STAR
(2) Unannot EA vs Unannot+Novel STAR

Look at: all possible, PB junctions only, Junction in genes with PB transcripts

*****/
%include "!MCLAB/event_analysis/sas_programs/analysis_junc_counts_star_vs_events_novel.sas";


/*


I think the true negative rate for EA junctions is inflated. I think I need filter
the unannotated catalog junctions. Try counting only those unannotated junctions
that have an APN>0 in at least one sample: keep only junctions that have at least
one read aligning to them
*/








