
ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Checks for BLAST results */

%macro blastCheck(xscripts);

data blast_results;
   set event.dtct_event_&xscripts._blast_results;
run;

/* Check : how many alignments have <80% identity? */

data flag_ident;
   set blast_results;
   if perc_identity < 90 then flag_perc_identity_lt_90=1;
   else flag_perc_identity_lt_90=0;
run;

proc freq data=flag_ident;
   tables flag_perc_identity_lt_90;
run; 

/* Check : how many alignments have a length <50 bp? */

data flag_length;
   set blast_results;
   if length < 50 then flag_length_lt_50=1;
   else flag_length_lt_50=0;
run;

proc freq data=flag_length;
   tables flag_length_lt_50;
run;

%mend;

%blastCheck(pacbio);
%blastCheck(refseq);


/* Pacbio BLAST: Split Pacbio ID into transcript_id, known or novel, match type */

data pb_blast_parse_pb;
  length pb_transcript_id $12.;
  length pb_status $8.;
  length pb_type $30.;
  set event.dtct_event_pacbio_blast_results;
  pb_transcript_id=scan(pacbio_id,1,"|");
  pb_status=scan(pacbio_id,2,"|");
  pb_type=scan(pacbio_id,3,"|");
  query_length=abs(query_stop-query_start)+1;
  ref_length=abs(ref_stop-ref_start)+1;
run;


/* Criteria for "hit":

   1. Percent identity > 95%
   2. No mismatches or gaps
   3. Length of alignment must be at least 95% the length of the full event sequence

*/

data event_len;
  set evspl.splicing_events_annot_refseq;
  keep event_id event_size transcript_id flag_junction_annotated flaG_intron_retention;
run;

data events_on;
   set event.splicing_on_apn_gt0;
   where flag_splicing_on=1;
   keep event_id;
run;


%macro bestHits(xscripts);

data blast_results;
   %if &xscripts.=pacbio %then %do;
   set pb_blast_parse_pb;
   %end; 
    %else %do;
   set event.dtct_event_&xscripts._blast_results;
   %end;
run;

proc sort data=blast_results;
  by event_id;
proc sort data=event_len;
  by event_id;
run;

data blast_results_w_len;
  merge blast_results (in=in1) event_len (in=in2);
  by event_id;
  if in1 and in2;
run;

data best_blast_hits;
   set blast_results_w_len;
   if mismatch > 0 then delete; *remove results with mismatches;
   if gapopen > 0 then delete; *remove results with gaps;
   if length ge 0.95*event_size then output;
run;

proc sort data=events_on;
   by event_id;
proc sort data=best_blast_hits;
   by event_id;
run;

data blast_hits_2_on;
  merge best_blast_hits (in=in1) events_on (in=in2);
  by event_id;
  if in1 then flag_blast_hit=1; else flag_blast_hit=0;
  if in2 then flag_splicing_on=1; else flag_splicing_on=0;
run;

*and just to check I have only detected events;
proc freq data=blast_hits_2_on;
   tables flag_blast_hit*flag_splicing_on;
run;

/* Make permenant -- I will cat all "best hits" per event together, flagging if they match to transcripts or not
   then compare whether this in agreement with their annotation */

data event.events_best_&xscripts._hit;
   set best_blast_hits;
run;

%mend;


%bestHits(pacbio);
%bestHits(refseq);
