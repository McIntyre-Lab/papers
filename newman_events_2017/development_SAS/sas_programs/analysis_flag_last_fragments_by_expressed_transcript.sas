ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* 3' variability check: check fragments of 3'-most exons matching to transcripts

For each transcript:
1. flag the fusion containing its last exon
2. flag the fragments of its last exon and number the order of fragments 5'->3' relative to orientation of transcript
3. flag last fragment
4. flag transcript if no detection in its last fusion
5. If last fusion detected, flag each fragment if not detected and cat together fragment ## detected and not detected

So, per transcript:

transcript
last fusion
last fragments
last fragment
last fusion dtct
# last fragments total
# last fragments dtct
cat fragments dtct (5>3)
cat fragments not dtct (5>3)

If all fragments are on, then no 3' variability
If some fragments off, then 3' variability

I also want to compare this to BLAST hits for PacBio later (need to think on how to do this)

I will probably later limit this to a subset of transcripts (e.g. ones with unique pieces)
but for now I will work on everything with at least one feature detected
*/

/* For each transcript, order the fragments 5'->3' and number */

/* drop short fragments */

data long_frags;
  set event.flagged_fragment_length;
  where flag_fragment_lt_min_bp=0;
  keep fragment_id;
run;

data frag2xs;
   set event.last_frag_exon_fus_by_xscript;
run;

proc sort data=long_frags;
  by fragment_id;
proc sort data=frag2xs;
  by fragment_id;
run;

data drop_short_frags;
   merge frag2xs (in=in1) long_frags (in=in2);
   by fragment_id;
   if in1 and in2;
run;

* split on strand;

data plus_xs minus_xs oops;
    set drop_short_frags;
    if strand="+" then output plus_xs;
    else if strand="-" then output minus_xs;
    else output oops;
run;

proc sort data=plus_xs;
  by transcript_id fragment_start fragment_end;
run;

proc sort data=minus_xs;
  by transcript_id descending fragment_end descending fragment_start;
run;

data plus_frag_num;
   set plus_xs;
   retain frag_num;
   by transcript_id;
   if first.transcript_id then frag_num=1;
   else frag_num=frag_num + 1;
   if last.transcript_id then flag_frag_last_in_region=1; else flag_frag_last_in_region=0;
run;

data minus_frag_num;
   set minus_xs;
   retain frag_num;
   by transcript_id;
   if first.transcript_id then frag_num=1;
   else frag_num=frag_num + 1;
   if last.transcript_id then flag_frag_last_in_region=1; else flag_frag_last_in_region=0;
run;

/* Flag fragment if last for transcript */

data plus_xs_frags_only;
   set plus_frag_num;
   where flag_frag_of_xs=1;
run;

data plus_flag_last_ex_frag;
   set plus_xs_frags_only;
   by transcript_id;
   if last.transcript_id then flag_frag_last_in_xscript=1;
   else flag_frag_last_in_xscript=0;
   keep transcript_id fusion_id fragment_id flag_frag_last_in_xscript;
run;

proc sort data=plus_frag_num;
   by transcript_id fusion_id fragment_id;
proc sort data=plus_flag_last_ex_frag;
   by transcript_id fusion_id fragment_id;
run;

data plus_frag_flag_last;
  merge plus_frag_num (in=in1) plus_flag_last_ex_frag;
   by transcript_id fusion_id fragment_id;
  if in1;
run;

data minus_xs_frags_only;
   set minus_frag_num;
   where flag_frag_of_xs=1;
run;

data minus_flag_last_ex_frag;
   set minus_xs_frags_only;
   by transcript_id;
   if last.transcript_id then flag_frag_last_in_xscript=1;
   else flag_frag_last_in_xscript=0;
   keep transcript_id fusion_id fragment_id flag_frag_last_in_xscript;
run;


proc sort data=minus_frag_num;
   by transcript_id fusion_id fragment_id;
proc sort data=minus_flag_last_ex_frag;
   by transcript_id fusion_id fragment_id;
run;

data minus_frag_flag_last;
  merge minus_frag_num (in=in1) minus_flag_last_ex_frag;
   by transcript_id fusion_id fragment_id;
  if in1;
run;

/* Stack together and make permenant */

data all_xs_flag_last_frag;
    set plus_frag_flag_last minus_frag_flag_last;
run;

data event.flag_last_frag_by_xscript;
   set all_xs_flag_last_frag;
   drop chrom start stop strand fragment_start fragment_end;
run;

