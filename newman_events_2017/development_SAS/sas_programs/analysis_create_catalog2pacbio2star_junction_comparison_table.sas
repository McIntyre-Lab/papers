
/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* I want to create a "master" dataset that relates catalog junctions with PacBio junctions and STAR junctions, including
   detection flags. I want to use this to figure out the false positive/negative rate with STAR/events relative to
   the PacBio junction

   Variables I will likely need:
   Junction coord		easiest to merge junctions on this
   flag_star_junction		is the junction detected with STAR?
   flag_catalog_junction	is the junction in the junction catalog?
   flag_pacbio_junction		is the junction seen in PacBio transcripts?
   flag_catalog_annotated_junc	is the junction considered "annotated"?
   flag_pacbio_novel_donor	does the PB junction have a novel donor site?
   flag_pacbio_novel_acceptor	does the PB junction have a novel acceptor site?
   flag_star_novel_donor	does the STAR junction have a novel donor site?
   flag_star_novel_acceptor	does the STAR junction have a novel acceptor site?
   STAR detection flags: 	detected in 1 sample (unique reads>0, >=5), 2 samples
   Events detection flags: 	detected in 1 sample (APN>0, >=5), 2 samples
*/

/* STAR junctions */
data star_junc;
  set event.npc_star_junctions;
  donor_stop=intron_start-1;
  acceptor_start=intron_stop;
  if strand=1 then strand_str="+";
  else strand_str="-";
  if NSC1 > 0 then flag_STAR_depth_NSC1_gt0=1; else flag_STAR_depth_NSC1_gt0=0;
  if NSC2 > 0 then flag_STAR_depth_NSC2_gt0=1; else flag_STAR_depth_NSC2_gt0=0;
  if  flag_STAR_depth_NSC1_gt0=1 and flag_STAR_depth_NSC2_gt0=1 then flag_STAR_depth_NSC_both_gt0=1;
  else flag_STAR_depth_NSC_both_gt=0;
  keep chr donor_stop acceptor_start strand_str
       flag_STAR_depth_NSC1_gt0 flag_STAR_depth_NSC2_gt0 flag_STAR_depth_NSC_both_gt0
       flag_depth_NSC1_ge5 flag_depth_NSC2_ge5 flag_NSC_all_ge5;
  rename strand_str=strand
         flag_depth_NSC1_ge5=flag_STAR_depth_NSC1_ge5
         flag_depth_NSC2_ge5=flag_STAR_depth_NSC2_ge5
         flag_NSC_all_ge5=flag_STAR_depth_NSC_both_ge5;
run;

/* PacBio junctions -- just the ones that are observed in at least one transcript 
   "unannotated" and border junctions here aren't useful ...
*/
data pb_junc;
   set evspl.splicing_events_annotations;
   where num_transcripts > 0;
   keep chr strand feature1_stop feature2_start;
   rename feature1_stop=donor_stop feature2_start=acceptor_start;
run;

/* Catalog junctions */
data cat_junc;
  set evspl.splicing_events_annot_refseq;
  if num_transcripts>0 then flag_junction_annotated=1;
  keep event_id chr strand feature1_stop feature2_start flag_intron_retention flag_junction_annotated;
  rename feature1_stop=donor_stop feature2_start=acceptor_start;
run;

/* Junction detection flags -- since I want these for NSC1 and NSC2 separately, I will need to
   reflag these */

data nsc1 nsc2;
   set event.mm10_refseq_splicing_counts;
   if sample_id="NSC1" then output nsc1;
   if sample_id="NSC2" then output nsc2;
   keep event_id apn;
run;

data flag_nsc1;
  set nsc1;
  if apn > 0 then flag_apn_NSC1_gt0=1; else flag_apn_NSC1_gt0=0;
  if apn ge 5 then flag_apn_NSC1_ge5=1; else flag_apn_NSC1_ge5=0;
  drop apn;
run;

data flag_nsc2;
  set nsc2;
  if apn > 0 then flag_apn_NSC2_gt0=1; else flag_apn_NSC2_gt0=0;
  if apn ge 5 then flag_apn_NSC2_ge5=1; else flag_apn_NSC2_ge5=0;
  drop apn;
run;

proc sort data=flag_nsc1;
  by event_id;
proc sort data=flag_nsc2;
  by event_id;
proc sort data=cat_junc;
  by event_id;
run;

data cat_junc_w_flags;
  merge cat_junc flag_nsc1 flag_nsc2;
  by event_id;
run;

/* PB junc -- I want to flag if the donor/acceptor site is seen in the RefSeq set
   This is for NNC junctions, as I want to track what parts are observable */

data cat_donors;
  set cat_junc_w_flags;
  where event_id ^? "intron";
  keep chr strand donor_stop;
run;

data cat_acceptors;
  set cat_junc_w_flags;
  where event_id ^? "intron";
  keep chr strand acceptor_start;
run;

proc sort data=cat_donors nodup;
  by chr strand donor_stop;
proc sort data=pb_junc nodup;
  by chr strand donor_stop acceptor_start;
run;

data pb_junc_flag_donor;
  merge pb_junc (in=in1) cat_donors (in=in2);
  by chr strand donor_stop;
  if in2 then flag_donor_in_catalog=1; else flag_donor_in_catalog=0;
  if in1 then output;
run;

proc sort data=pb_junc_flag_donor;
   by chr strand acceptor_start;
proc sort data=cat_acceptors nodup;
   by chr strand acceptor_start;
run;

data pb_junc_flag_acceptor;
  merge pb_junc_flag_donor (in=in1) cat_acceptors (in=in2);
  by chr strand acceptor_start;
  if in2 then flag_acceptor_in_catalog=1; else flag_acceptor_in_catalog=0;
  if in1 then output;
run;

/* Merge PB junctions and catalog junctions */

proc sort data=cat_junc_w_flags;
   by chr strand donor_stop acceptor_start;
proc sort data=pb_junc_flag_acceptor;
   by chr strand donor_stop acceptor_start;
run;

data cat_junc_w_pb;
  length chr $20.;
  format chr $20.;
  merge cat_junc_w_flags (in=in1) pb_junc_flag_acceptor (in=in2);
  by chr strand donor_stop acceptor_start;
  if in1 then flag_in_catalog=1; else flag_in_catalog=0;
  if in2 then flag_in_pacbio=1; else flag_in_pacbio=0;
run;

proc sort data=star_junc;
   by chr strand donor_stop acceptor_start;
run;

data junc_compare;
  merge cat_junc_w_pb (in=in1) star_junc (in=in2);
  by chr strand donor_stop acceptor_start ;
  if in2 then flag_in_star=1; else flag_in_star=0;
  if not in1 then do;
      flag_in_pacbio=0;
      flag_in_catalog=0;
      end;
run;

/* Make permenant */

data event.catalog_pacbio_star_junctions;
   set junc_compare;
run;


