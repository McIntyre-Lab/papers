
ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* PacBio BLAST -- For each event, I am going to cat together the PacBio transcripts they map to.
   Flag if at least one PB transcript is "Known" and at least one is "Novel"
   Then flag if at least transcript is NIC, NNC, etc. */

proc freq data=event.events_best_pacbio_hit;
   tables pb_status pb_type;
run;

/*
PB types:
full-splice_match		all junctions of reference transcript match
incomplete-splice_match		consective but not all junctions of ref transcript match
novel_in_catalog		novel transcript, existing donors/acceptors
novel_not_in_catalog		novel transcript, novel donors/acceptors
intergenic			transcript placed outside the boundaries of an known gene
genic_intron			transcript placed completely within an intronic region of known gene
genic				monoexonic transcript, overlapping existing exons and introns
fusion				fusion transcript
antisense			transcript produced from opposite strand of annotated gene
*/

data pb_blast;
  set event.events_best_pacbio_hit;

  /* Flag known/novel PB transcripts */
  if pb_status="Known" then flag_pb_transcript_known=1; else flag_pb_transcript_known=0;
  if pb_status="Novel" then flag_pb_transcript_novel=1; else flag_pb_transcript_novel=0;

  /* Flag PB transcript type */
  if pb_type="antisense" then flag_pb_antisense=1; else flag_pb_antisense=0;
  if pb_type="full-splice_match" then flag_pb_fsm=1; else flag_pb_fsm=0;
  if pb_type="fusion" then flag_pb_fusion=1; else flag_pb_fusion=0;
  if pb_type="genic" then flag_pb_genic=1; else flag_pb_genic=0;
  if pb_type="genic_intron" then flag_pb_genic_intron=1; else flag_pb_genic_intron=0;
  if pb_type="incomplete-splice_match" then flag_pb_ism=1; else flag_pb_ism=0;
  if pb_type="intergenic" then flag_pb_intergenic=1; else flag_pb_intergenic=0;
  if pb_type="novel_in_catalog" then flag_pb_nic=1; else flag_pb_nic=0;
  if pb_type="novel_not_in_catalog" then flag_pb_nnc=1; else flag_pb_nnc=0;

  if transcript_id ne '' then flag_junction_annotated=1;
  drop pacbio_id perc_identity length mismatch gapopen query_start query_stop
       ref_start ref_stop evalue bitscore query_length ref_length event_size;
run;

/* Sum flags */

proc sort data=pb_blast nodup;
   by event_id pb_transcript_id;
run;

proc means data=pb_blast noprint;
   by event_id;
   var flag_pb_transcript_known flag_pb_transcript_novel
       flag_pb_antisense flag_pb_fsm flag_pb_fusion flag_pb_genic flag_pb_genic_intron
       flag_pb_ism flag_pb_intergenic flag_pb_nic flag_pb_nnc
       flag_junction_annotated flag_intron_retention;
   output out=pb_blast_sum_hits
          sum(flag_pb_transcript_known)=num_pb_known
          sum(flag_pb_transcript_novel)=num_pb_novel
          sum(flag_pb_antisense)=num_pb_novel_antisense
          sum(flag_pb_fsm)=num_pb_full_match
          sum(flag_pb_fusion)=num_pb_novel_fusion
          sum(flag_pb_genic)=num_pb_novel_genic
          sum(flag_pb_genic_intron)=num_pb_novel_genic_intron
          sum(flag_pb_ism)=num_pb_incomplete_match
          sum(flag_pb_intergenic)=num_pb_novel_intergenic
          sum(flag_pb_nic)=num_pb_novel_in_catalog
          sum(flag_pb_nnc)=num_pb_novel_not_catalog
          max(flag_pb_transcript_known)=flag_pb_known
          max(flag_pb_transcript_novel)=flag_pb_novel
          max(flag_pb_antisense)=flag_pb_novel_antisense
          max(flag_pb_fsm)=flag_pb_full_match
          max(flag_pb_fusion)=flag_pb_novel_fusion
          max(flag_pb_genic)=flag_pb_novel_genic
          max(flag_pb_genic_intron)=flag_pb_novel_genic_intron
          max(flag_pb_ism)=flag_pb_incomplete_match
          max(flag_pb_intergenic)=flag_pb_novel_intergenic
          max(flag_pb_nic)=flag_pb_novel_in_catalog
          max(flag_pb_nnc)=flag_pb_novel_not_catalog
          max(flag_junction_annotated)=flag_junction_annotated
          max(flag_intron_retention)=flag_intron_retention;
run;

data pb_xs;
  set pb_blast;
  keep event_id pb_transcript_id;
run;

proc sort data=pb_xs nodup;
  by event_id pb_transcript_id;
run;
proc freq noprint data=pb_xs;
   tables event_id / out=xs_count;
proc sort data=xs_count;
   by descending count;
proc print data=xs_count (obs=1);
run; *31 xs;

data cat_xs;
   array xs[31] $15.;
   retain xs1-xs31;
   set pb_xs;
   by event_id;
   if first.event_id then do;
      call missing(of xs1-xs31);
      records = 0;
   end;
   records + 1;
   xs[records]=pb_transcript_id;
   if last.event_id then output;
run;

data cat_xs2;
  set cat_xs;
  length pb_transcript_id2 $400.;
  rename records=num_pb_transcripts;
  pb_transcript_id2=catx("|", OF xs1-xs31);
  drop xs1-xs31 pb_transcript_id;
  rename pb_transcript_id2=pb_transcript_id;
run;

proc sort data=cat_xs2;
   by event_id;
proc sort data=pb_blast_sum_hits;
   by event_id;
run;

data pb_blast_w_annot;
   merge cat_xs2 (in=in1) pb_blast_sum_hits (in=in2);
   by event_id;
   if in1 and in2;
   drop _TYPE_ _FREQ_;
run;

/* Make permenant */


data event.pacbio_blast_hits_w_annot;
  set pb_blast_w_annot;
run;

/* Count flags */

proc freq data=pb_blast_w_annot noprint;
    tables flag_junction_annotated*flag_intron_retention*flag_pb_known*flag_pb_novel / out=pb_checks;
run;

proc print data=pb_checks;
run;

/*

 flag_junction_    flag_intron_    flag_pb_    flag_pb_
    annotated        retention       known       novel     COUNT    PERCENT

        0                0             0           1         214     0.3358
        0                0             1           0          18     0.0282
        0                0             1           1           3     0.0047
        0                1             0           1         997     1.5647
        0                1             1           0         263     0.4127
        0                1             1           1          95     0.1491
        1                0             0           1        3685     5.7831
        1                0             1           0       36803    57.7574
        1                0             1           1       21642    33.9642

Some annotated junctions going to novel PB transcript (this is okay)
Unannotated junctions and IR eventsgoing to known transcripts (??)
*/



