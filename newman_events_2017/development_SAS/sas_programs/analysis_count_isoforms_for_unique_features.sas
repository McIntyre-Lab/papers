/* Summarizing some event analysis stuff:

(1) Count the number of detected events and fragments:
	(a) by type (fusion, fragment, junction, unannotated junc, etc.)
	(b) by commonality type (unique, common, constitutive)

(2) For unique features detected, how many isoforms?
	(a) fragments only
	(b) fusions only
	(c) events only
	(d) fragments + events
	(e) fusions + events
	(f) overlap between (a) and (b)
	(g) overlap between (a) and (c)
	(h) overlap between (b) and (c)
	(i) overlap between (d) and (e)

(3) Calculate the distribution of transcripts by number of detected unique pieces
	(a) #uniq * transcript (features, fragments, events, fusions)
	(b) %uniq * transcript (features, fragments, events, fusions)

(4) Bin transcripts by whether they have a unique pieces, common or constitutive
	e.g. # xs with at least one uniq detected
	     # xs with no uniq, but have common detected
             # xs with no uniq, no common, only constitutive
	(a) fusions + events
	(b) fragments + events

(ignore multigene for now)
*/

libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For unique features detected, how many isoforms?
	(a) fragments only
	(b) fusions only
	(c) events only
	(d) fragments + events
	(e) fusions + events
	(f) overlap between (a) and (b)
	(g) overlap between (a) and (c)
	(h) overlap between (b) and (c)
	(i) overlap between (d) and (e) */

data unique_features;
  set event.dtct_features_by_class;
  where flag_feature_unique=1;
  keep feature_id feature_type;
run;

data frag2xs;
   set mm10.mm10_exon_fragment_flagged;
   where num_xscripts_per_fragment=1;
   keep fragment_id transcript_id;
   rename fragment_id=feature_id;
run;

data fusion2xs;
   set mm10.mm10_si_fusions_unique_flagged;
   where num_xscripts_per_fusion=1;
   keep fusion_id transcript_id;
   rename fusion_id=feature_id;
run;

data event2xs;
   set evspl.splicing_events_annot_refseq;
   where num_transcripts = 1;
   keep event_id transcript_id;
   rename event_id=feature_id;
run;

data feature2xs;
  set event2xs fusion2xs frag2xs;
run;

proc sort data=feature2xs;
   by feature_id;
proc sort data=unique_features;
   by feature_id;
run;

data uniq_feat2xs;
  merge unique_features (in=in1) feature2xs (in=in2);
  by feature_id;
  if in1 and in2;
run;

data xs_w_uniq_frag  xs_w_uniq_fus  xs_w_uniq_splice oops;
   set uniq_feat2xs;
   if feature_type="fragment" then output xs_w_uniq_frag;
   else if feature_type="fusion" then output xs_w_uniq_fus;
   else if feature_type="splicing" then output xs_w_uniq_splice;
   else output oops;
   keep transcript_id;
run;

proc sort data=xs_w_uniq_frag nodup;
   by transcript_id;
proc sort data=xs_w_uniq_fus nodup;
   by transcript_id;
proc sort data=xs_w_uniq_splice nodup;
   by transcript_id;
run;

data event.flag_xscript_uniq_features;
  merge xs_w_uniq_frag (in=in1) xs_w_uniq_fus (in=in2) xs_w_uniq_splice (in=in3);
  by transcript_id;
  if in1 then flag_xs_has_uniq_fragment=1; else flag_xs_has_uniq_fragment=0;
  if in2 then flag_xs_has_uniq_fusion=1; else flag_xs_has_uniq_fusion=0;
  if in3 then flag_xs_has_uniq_event=1; else flag_xs_has_uniq_event=0;
  if in1 and in3 then flag_xs_has_uniq_fragment_event=1; else flag_xs_has_uniq_fragment_event=0;
  if in2 and in3 then flag_xs_has_uniq_fusion_event=1; else flag_xs_has_uniq_fusion_event=0;
run;

proc freq data=event.flag_xscript_uniq_features;
   tables flag_xs_has_uniq_fragment
          flag_xs_has_uniq_fusion
          flag_xs_has_uniq_event
          flag_xs_has_uniq_fragment*flag_xs_has_uniq_event
          flag_xs_has_uniq_fusion*flag_xs_has_uniq_event
          flag_xs_has_uniq_fragment_event*flag_xs_has_uniq_fusion_event;
run;

/*

 flag_xs_has_
uniq_fragment    Frequency
---------------------------
            0        5574
            1       24056


 flag_xs_has_
  uniq_fusion    Frequency
---------------------------
            0       15603
            1       14027


 flag_xs_has_
   uniq_event    Frequency
---------------------------
            0       15252
            1       14378

xs with unique:
fragment:	24056
fusion:		14027
splicing:	14378

Overlap between uniq_fragments and uniq_events:

 flag_xs_has_uniq_fragment
           flag_xs_has_uniq_event

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |    108 |   5466 |   5574
          |   0.36 |  18.45 |  18.81
          |   1.94 |  98.06 |
          |   0.71 |  38.02 |
 ---------+--------+--------+
        1 |  15144 |   8912 |  24056
          |  51.11 |  30.08 |  81.19
          |  62.95 |  37.05 |
          |  99.29 |  61.98 |
 ---------+--------+--------+
 Total       15252    14378    29630
             51.47    48.53   100.00

xs with uniq fragment only: 15144
xs with uniq event only: 5466
xs with uniq fragment and event: 8912

Overlap between uniq_fusions and uniq_events:

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |   8244 |   7359 |  15603
          |  27.82 |  24.84 |  52.66
          |  52.84 |  47.16 |
          |  54.05 |  51.18 |
 ---------+--------+--------+
        1 |   7008 |   7019 |  14027
          |  23.65 |  23.69 |  47.34
          |  49.96 |  50.04 |
          |  45.95 |  48.82 |
 ---------+--------+--------+
 Total       15252    14378    29630
             51.47    48.53   100.00

xs with uniq fusion only: 7008
xs with uniq event only: 7359
xs with uniq fusion and event: 7019

*/
