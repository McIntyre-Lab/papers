
ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* RefSeq BLAST -- For each event, I am going to cat together the set of transcripts it maps to,
   and cross-check with event type (annotated/unannotated junction, IR event).
   Then, if an event is already annotated, check that it maps to the same set of transcripts */


data event_annot;
  set event.events_best_refseq_hit;
  if transcript_id ne '' then flag_junction_annotated=1;
  keep event_id transcript_id flag_junction_annotated flag_intron_retention;
run;

data event_hit;
  set event.events_best_refseq_hit;
  keep event_id refseq_id;
run;

data event2xs;
  length refseq_id $15.;
  set event_annot;
  if flag_junction_annotated=1 or transcript_id ne "";
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "") ;
     refseq_id=scan(transcript_id,i,"|");
     output;
     end;
  keep event_id refseq_id;
run;

proc sort data=event_hit nodup;
  by event_id refseq_id;
proc sort data=event2xs nodup;
  by event_id refseq_id;
run;

data both_xs  annot_xs blast_xs;
  merge event2xs (in=in1) event_hit (in=in2);
  by event_id refseq_id;
  if in1 and in2 then output both_xs; *same number going in (568674), so we get all the annotations back;
  else if in1 then output annot_xs; *0 obs!;
  else output blast_xs; *3347 obs. So some events are mapping to additional transcripts?;
run;


/* Cat together the set of BLAST-only transcripts  */

proc sort data=blast_xs;
   by event_id refseq_id;
proc freq noprint data=blast_xs;
   tables event_id / out=xs_count;
proc sort data=xs_count;
   by descending count;
proc print data=xs_count (obs=1);
run; *120 xs??;

data cat_xs;
   array xs[120] $15.;
   retain xs1-xs120;
   set blast_xs;
   by event_id;
   if first.event_id then do;
      call missing(of xs1-xs120);
      records = 0;
   end;
   records + 1;
   xs[records]=refseq_id;
   if last.event_id then output;
run;

data cat_xs2;
  set cat_xs;
  length blast_only_transcript_id $1500.;
  rename records=num_blast_only_transcripts;
  blast_only_transcript_id=catx("|", OF xs1-xs120);
  drop xs1-xs120 refseq_id;
run;

data both_event_only;
  set both_xs;
  keep event_id;
run;


proc sort data=cat_xs2;
  by event_id;
proc sort data=both_event_only nodup;
  by event_id;
proc sort data=event_annot;
  by event_id;
run;

data event_annot_w_blast_hits;
  merge event_annot (in=in1) both_event_only (in=in2) cat_xs2 (in=in3);
  by event_id;
  if in2 or in3 then flag_blast_hit=1; else flag_blast_hit=0;
  if in2 then flag_annotated_blast_hit=1; else flag_annotated_blast_hit=0;
  if in3 then flag_unannotated_blast_hit=1; else flag_unannotated_blast_hit=0;
run; 


/* Make permenant */

data event.refseq_blast_hits_w_annot;
  set event_annot_w_blast_hits;
run;


proc freq data=event_annot_w_blast_hits noprint;
  tables flag_junction_annotated*flag_intron_retention*
         flag_blast_hit*flag_annotated_blast_hit*flag_unannotated_blast_hit / out = flag_check;
run;

proc print data=flag_check;
run;


/*

                               flag_
flag_junction_  flag_intron_  blast_  flag_annotated_  flag_unannotated_
   annotated      retention     hit      blast_hit         blast_hit       COUNT  PERCENT

       0              0          1           0                 1             121   0.0211
       0              1          1           0                 1            1681   0.2936
       1              0          1           0                 1              32   0.0056
       1              0          1           1                 0          566624  98.9489
       1              0          1           1                 1            4185   0.7308


Vast majority of hits are annotated junctions to their annotated transcripts, so this is good.

There are:
121 unannotated junctions map to transcripts
1681 IR events maps to transcripts
32 annotated junctions that only map to transcripts they aren't annotated to
4185 annotated junctions that map to their associated transcripts as well as additional transcripts

Hypothesis: annotation problem??
*/

  

