ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For the set of RefSeq transcripts with all features detected,
   how many have a PacBio match?
   Look at 100% on, 75% on, 50% on, 25% on */

data xs_100_on xs_75_on xs_50_on xs_25_on;
   set event.xscripts_w_unique_by_bin;
   if perc_features_dtct = 1 then output xs_100_on;
   if perc_features_dtct ge 0.75 then output xs_75_on;
   if perc_features_dtct ge 0.5 then output xs_50_on;
   if perc_features_dtct ge 0.25 then output xs_25_on;
   keep transcript_id;
run;

*19338 with all on;
*48913 with 75%+ on;
*65412 with 50%+ on;
*81165 with 25%+ on;

data pb2refseq;
   set splice_match;
   where transcript_id ^? "ENSMUST";
   length splice_match_id $18.;
   splice_match_id=scan(transcript_id,1,".");
   keep splice_match_id match_type;
   rename splice_match_id=transcript_id;
run; *8722;

proc freq data=pb2refseq;
   tables match_type;
run;


proc sort data=pb2refseq nodup;
  by transcript_id;
proc sort data=xs_100_on;
  by transcript_id;
proc sort data=xs_75_on;
  by transcript_id;
proc sort data=xs_50_on;
  by transcript_id;
proc sort data=xs_25_on;
  by transcript_id;
run;


data xs_100_on_pb;
  merge xs_100_on (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;

proc freq data=xs_100_on_pb;
   where match_type ne "incomplete_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

* 13290 RefSeq without PB, 6048 in both, 2167 PB no Refseq 100%;

* 13290 Refseq without PB, 5669 in both, 1339 full-spliced matched PB no Refseq (81%);

proc freq data=xs_100_on_pb;
   where match_type ne "full_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

* 13290 Refseq without PB, 615 in both, 869 incomplete-spliced matched PB no Refseq (41%);


data xs_75_on_pb;
  merge xs_75_on (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;

proc freq data=xs_75_on_pb;
   where match_type ne "incomplete_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_75_on_pb;
   where match_type ne "full_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;
* 41100 RefSeq without PB, 7813 in both, 402 PB no Refseq 100%;

* 41100 Refseq without PB, 6715 in both, 293 full-spliced matched PB no Refseq (96%);
* 41100 Refseq without PB, 1374 in both, 110 incomplete-spliced matched PB no Refseq (93%);

data xs_50_on_pb;
  merge xs_50_on (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;


proc freq data=xs_50_on_pb;
   where match_type ne "incomplete_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_50_on_pb;
   where match_type ne "full_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;
* 57331 RefSeq without PB, 8081 in both, 134 PB no Refseq 100%;

* 57331 Refseq without PB, 6888 in both, 120 full-spliced matched PB no Refseq (98%);
* 57331 Refseq without PB, 1470 in both, 14 incomplete-spliced matched PB no Refseq (99%);


data xs_25_on_pb;
  merge xs_25_on (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;

proc freq data=xs_25_on_pb;
   where match_type ne "incomplete_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_25_on_pb;
   where match_type ne "full_splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;
* 73050 RefSeq without PB, 8115 in both, 100 PB no Refseq 100%;

* 73050 Refseq without PB, 6915 in both, 93 full-spliced matched PB no Refseq (99%);
* 73050 Refseq without PB, 1477 in both, 7 incomplete-spliced matched PB no Refseq (100%);


/* For the set of 93 full-spliced match PB transcripts that have no RefSeq transcript with at least 25% of total
   features detected, check to see if these contain the 83 PB transcripts associated with non-expressed transcripts */


data pb_check;
   set xs_25_on_pb;
   where flag_xs_has_pb=1 and flag_xs_dtct=0;
   keep transcript_id;
run;

data pb_refseq_off;
   set event.refseq_xscripts_with_pacbio;
   where flag_xscript_has_pacbio=1 and bin_xscript_perc_uniq_dtct="not_exp";
run;


proc sort data=pb_check;
   by transcript_id;
proc sort data=pb_refseq_off;
   by transcript_id;
run;

data pb_check2;
  merge pb_check (in=in1) pb_refseq_off (in=in2);
  by transcript_id;
  if in1 then flag_xscript_low_dtct=1; else flag_xscript_low_dtct=0;
  if in2 then flag_xscript_dropped=1; else flag_xscript_dropped=0;
run;

proc freq data=pb_check2;
   tables flag_xscript_low_dtct*flag_xscript_dropped;
run;

data bins;
set event.xscripts_w_unique_by_bin;
run;

proc sort data=bins;
   by transcript_id;
proc sort data=pb_refseq_off;
   by transcript_id;
run;

data check;
  merge bins (in=in1) pb_refseq_off (in=in2);
   by transcript_id;
  if in1 and in2;
run;


*the 83 PB transcripts are from genes we eliminated due to no expression of features;
