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

*14734 with all on;
*34622 with 75%+ on;
*45883 with 50%+ on;
*57572 with 25%+ on;

data pb2refseq;
   set event.pacbio2refseq_id_nomulti;
   if transcript_id="XM_011239942" then match_type="incomplete-splice_match";
   keep transcript_id match_type;
run;

proc sort data=pb2refseq nodup;
   by transcript_id;
run;

proc freq data=pb2refseq;
   tables match_type;
run;

/*

                                                     Cumulative
 match_type                 Frequency     Percent     Frequency
 --------------------------------------------------------------
 full-splice_match              5412       82.17          5412        82.17
 incomplete-splice_match        1174       17.83          6586       100.00



*/

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
   where match_type ne "incomplete-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

* 11433 RefSeq without PB, 3468 in both, 3118 PB no Refseq 100%;

* 11433 Refseq without PB, 3104 in both, 2308 full-spliced matched PB no Refseq (81%);

proc freq data=xs_100_on_pb;
   where match_type ne "full-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

* 11433 Refseq without PB, 364 in both, 810 incomplete-spliced matched PB no Refseq (41%);


data xs_75_on_pb;
  merge xs_75_on (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;


proc freq data=xs_75_on_pb;
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_75_on_pb;
   where match_type ne "incomplete-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_75_on_pb;
   where match_type ne "full-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;
* 30314 RefSeq without PB, 4504 in both, 2082 PB no Refseq 100%;

* 30314 Refseq without PB, 3757 in both, 1655 full-spliced matched PB no Refseq (96%);
* 30314 Refseq without PB, 747 in both, 427 incomplete-spliced matched PB no Refseq (93%);

data xs_50_on_pb;
  merge xs_50_on (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;


proc freq data=xs_50_on_pb;
      tables flag_xs_dtct*flag_xs_has_pb ;
run;


proc freq data=xs_50_on_pb;
   where match_type ne "incomplete-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_50_on_pb;
   where match_type ne "full-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;
* 41453 RefSeq without PB, 4626 in both, 1960 PB no Refseq 100%;

* 41453 Refseq without PB, 3827 in both, 1585 full-spliced matched PB no Refseq (98%);
* 41453 Refseq without PB, 799 in both, 375 incomplete-spliced matched PB no Refseq (99%);


data xs_25_on_pb;
  merge xs_25_on (in=in1) pb2refseq (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;


proc freq data=xs_25_on_pb;
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_25_on_pb;
   where match_type ne "incomplete-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;

proc freq data=xs_25_on_pb;
   where match_type ne "full-splice_match";
   tables flag_xs_dtct*flag_xs_has_pb ;
run;
* 53118 RefSeq without PB, 4650 in both, 1936 PB no Refseq 100%;
* 53118 Refseq without PB, 3848 in both, 1564 full-spliced matched PB no Refseq (99%);
* 53118 Refseq without PB, 802 in both, 372 incomplete-spliced matched PB no Refseq (100%);


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
