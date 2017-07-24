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


     data WORK.SPLICE_MATCH    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile '!MCLAB/event_analysis/references/pacbio_isoforms_list_for_import.txt'
 delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat pacbio_id $10. ; informat source $6.;    informat feature_type $10. ;
        informat start best32. ;  informat end best32.;   informat frame best32. ;
        informat strand $1. ;     informat score best32.; informat transcript_id $18. ;
        informat match_type $23.; informat note $28. ;
        format pacbio_id $10. ;   format source $6. ;     format feature_type $10. ;
        format start best12. ;    format end best12. ;    format frame best12. ;
        format strand $1. ;       format score best12. ;  format transcript_id $18. ;
        format match_type $23. ;  format note $28. ;
        input pacbio_id $ source $ feature_type $ start end
              frame strand $ score transcript_id $ match_type $ note $
              ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;

data pb2refseq;
   set splice_match;
   where transcript_id ^? "ENSMUST";
   length splice_match_id $18.;
   splice_match_id=scan(transcript_id,1,".");
   keep pacbio_id splice_match_id match_type;
   rename splice_match_id=transcript_id;
run; *8722;

data pb2keep;
  set event.pacbio2refseq_id_nomulti;
  keep pacbio_id;
run;

proc sort data=pb2keep nodup;
  by pacbio_id;
proc sort data=pb2refseq nodup;
  by pacbio_id transcript_id;
run;

data pb2refseq2;
  merge pb2refseq (in=in1) pb2keep (in=in2);
  by pacbio_id;
  if in1 and in2;
  drop pacbio_id;
run;

proc sort data=pb2refseq2 nodup;
   by transcript_id;
run;

proc freq data=pb2refseq2;
   tables match_type;
run;

/*
                                                     Cumulative    Cumulative
 match_type                 Frequency     Percent     Frequency      Percent
 ----------------------------------------------------------------------------
 full_splice_match              5078       83.20          5078        83.20
 incomplete_splice_match        1025       16.80          6103       100.00


*/


proc freq data=pb2refseq2 noprint;
  tables transcript_id / out=xs_count;
proc sort data=xs_count;
  by descending count;
proc print data=xs_count(obs=1);
run; *3 matches max;

data pb_full pb_partial;
   set pb2refseq2;
   if match_type="full_splice_match" then output pb_full;
   else output pb_partial;
   keep transcript_id;
run;

proc sort data=pb_full nodup;
   by transcript_id;
proc sort data=pb_partial nodup;
   by transcript_id;
run;

data pb2refseq3;
  merge pb_full (in=in1) pb_partial (in=in2);
  by transcript_id;
  if in1 then flag_full_splice_match=1;
  else flag_full_splice_match=0;
run;

proc freq data=pb2refseq3;
   tables flag_full_splice_match;
run;

/*
       flag_full_                             Cumulative    Cumulative
     splice_match    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0         802       13.64           802        13.64
                1        5078       86.36          5880       100.00


*/


proc sort data=pb2refseq3 nodup;
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
  merge xs_100_on (in=in1) pb2refseq3 (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;

proc freq data=xs_100_on_pb noprint;
   tables flag_xs_dtct*flag_xs_has_pb*flag_full_splice_match / out=xs_100_count ;
proc print data=xs_100_count;
run;

/*
                               flag_full_
       flag_xs_    flag_xs_      splice_
Obs      dtct       has_pb        match      COUNT    PERCENT

 1         0           1            0          539     9.1667
 2         0           1            1          803    13.6565
 3         1           0            .        10196      .
 4         1           1            0          263     4.4728
 5         1           1            1         4275    72.7041

% FSM with hit: 4275/(803+4275) = 84.2%
% PSM with hit: 263/(263+539) = 32.8%

*/

data xs_75_on_pb;
  merge xs_75_on (in=in1) pb2refseq3 (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;

proc freq data=xs_75_on_pb noprint;
   tables flag_xs_dtct*flag_xs_has_pb*flag_full_splice_match / out=xs_75_count ;
proc print data=xs_75_count;
run;

/*
                         flag_full_
 flag_xs_    flag_xs_      splice_
   dtct       has_pb        match      COUNT

     0           1            0           68
     0           1            1          126
     1           0            .        28936
     1           1            0          734
     1           1            1         4952

% FSM with hit: 4952 / (4952+126) = 97.5%
% PSM with hit: 734 / (734+68) = 91.5%

*/



data xs_50_on_pb;
  merge xs_50_on (in=in1) pb2refseq3 (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;

proc freq data=xs_50_on_pb noprint;
   tables flag_xs_dtct*flag_xs_has_pb*flag_full_splice_match / out=xs_50_count ;
proc print data=xs_50_count;
run;

/*
                          flag_full_
  flag_xs_    flag_xs_      splice_
    dtct       has_pb        match      COUNT

      0           1            0            6
      0           1            1           44
      1           0            .        40053
      1           1            0          796
      1           1            1         5034

% FSM with hit: 5034/(5034+44) = 99.1%
% PSM with hit: 796/(796+6) = 99.3%

*/



data xs_25_on_pb;
  merge xs_25_on (in=in1) pb2refseq3 (in=in2);
  by transcript_id;
  if in1 then flag_xs_dtct=1; else flag_xs_dtct=0;
  if in2 then flag_xs_has_pb=1; else flag_xs_has_pb=0;
run;

proc freq data=xs_25_on_pb noprint;
   tables flag_xs_dtct*flag_xs_has_pb*flag_full_splice_match / out=xs_25_count ;
proc print data=xs_25_count;
run;

/*
                         flag_full_
 flag_xs_    flag_xs_      splice_
   dtct       has_pb        match      COUNT

     0           1            0            2
     0           1            1           22
     1           0            .        51716
     1           1            0          800
     1           1            1         5056

% FSM with hit: 5056/(5056+22) = 99.6%
% PSM with hit: 800/(800+2) = 99.8%

*/



/* For the set of 24 full-spliced match PB transcripts that have no RefSeq transcript with at least 25% of total
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
