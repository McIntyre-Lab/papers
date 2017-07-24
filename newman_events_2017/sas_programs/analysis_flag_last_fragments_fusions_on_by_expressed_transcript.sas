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

I also want to compare this to BLAST hits for PacBio later (need to think on how to do this)

I will probably later limit this to a subset of transcripts (e.g. ones with unique pieces)
but for now I will work on everything with at least one feature detected
*/

/* Merge in "on" flags for fusions and fragments */

data frag_on;
    set event.fragments_on_apn_gt0;
run;

data fus_on;
   set event.fusions_on_apn_gt0;
run;

data frag2xs;
  set event.flag_last_frag_by_xscript;
run;

proc sort data=frag_on;
   by fragment_id;
proc sort data=frag2xs;
   by fragment_id;
run;

data frag2xs_frag_on;
  merge frag2xs (in=in1) frag_on (in=in2);
  by fragment_id;
  if in1 and in2;
run;

proc sort data=frag2xs_frag_on;
   by fusion_id;
proc sort data=fus_on;
   by fusion_id;
run;

data frag2xs_frag_fus_on;
   merge frag2xs_frag_on (in=in1) fus_on (in=in2);
   by fusion_id;
   if in1 and in2;
run;

/* Count: how many transcripts have their 3'-most fusion not detected? */

data fus2xs;
  set frag2xs_frag_fus_on;
  keep transcript_id fusion_id flag_fusion_on;
run;

proc sort data=fus2xs nodup;
  by transcript_id fusion_id;
run;

proc freq data=fus2xs;
  tables flag_fusion_on;
run;

/*
                                             Cumulative    Cumulative
  flag_fusion_on    Frequency     Percent     Frequency      Percent
  -------------------------------------------------------------------
               0        9721       13.22          9721        13.22
               1       63807       86.78         73528       100.00


9721 transcripts where the 3'-most fusion is not detected.
This could be 3' variation, or it could be that the transcript isn't expressed
Exclude these from analysis

Focus only on the 63807 transcripts instead.
*/

/* make permenant */

data event.xscripts_last_fus_frag_onflags;
   set frag2xs_frag_fus_on;
run;


