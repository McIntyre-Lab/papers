ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* PACBIO BLAST PROCESSING:
Flag if events/fragments align to PB transcripts
Flag if events/fragments align to PB transcripts corresponding to at least one of their originating transcripts
Flag if known, novel PB
*/

/* EVENTS */

*non-redundant events;
data nr_events;
  set event.flag_redundant_event_rfsq_blast;
  keep event_id flag_event_redundant;
run;

*BLAST hits;

data blast_hits;
  length pacbio_id $10.;
  set event.events_best_pacbio_hit;
  keep pb_transcript_id pb_status pb_type event_id transcript_id flag_junction_annotated flag_intron_retention;
  rename pb_transcript_id=pacbio_id transcript_id=transcripts_cat;
run;

* PB-to-Refseq;

data pb2refseq;
  set event.pacbio2refseq_id;
  keep pacbio_id transcript_id;
run;

proc sort data=pb2refseq;
   by pacbio_id;
proc sort data=blast_hits;
   by pacbio_id;
run;

data blast_hits_pb2rs;
  merge blast_hits (in=in1) pb2refseq (in=in2);
  by pacbio_id;
  if in1 and in2 then do;
     flag_pb_has_refseq=1;
     output;
     end;
  else if in1 then do;
     flag_pb_has_refseq=0;
     output;
     end;
run;

proc sort data=blast_hits_pb2rs;
  by event_id;
proc sort data=nr_events nodup;
  by event_id;
run;

data blast_hits_pb2rs_nr;
  merge blast_hits_pb2rs (in=in1) nr_events (in=in2);
  by event_id;
  if in1 and in2 then output;
  else if in1 then do;
     flag_event_redundant=0;
     output;
     end;
run;

*check;

proc freq data=blast_hits_pb2rs_nr;
  tables flag_pb_has_refseq*flag_event_redundant;
run;

/*
  flag_pb_has_refseq
            flag_event_redundant

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  47670 |    128 |  47798
           |  35.99 |   0.10 |  36.09
           |  99.73 |   0.27 |
           |  36.05 |  61.84 |
  ---------+--------+--------+
         1 |  84559 |     79 |  84638
           |  63.85 |   0.06 |  63.91
           |  99.91 |   0.09 |
           |  63.95 |  38.16 |
  ---------+--------+--------+
  Total      132229      207   132436
              99.84     0.16   100.00
*/

/* Flag if PB transcript is in RefSeq list */

data flag_pb_in_refseq;
   set blast_hits_pb2rs_nr;
   if flag_pb_has_refseq=1 then do;
   if count(transcripts_cat,strip(transcript_id)) > 0 then flag_pb_in_refseq_list=1;
   else flag_pb_in_refseq_list=0;
   end;
run;

/* See if event aligns to any RefSeq transcript */

proc sort data=flag_pb_in_refseq;
   by event_id;
run;

proc means data=flag_pb_in_refseq noprint;
   by event_id;
   var flag_junction_annotated flag_intron_retention flag_pb_has_refseq flag_event_redundant
       flag_pb_in_refseq_list;
   output out=flag_pb_by_event max=;
run;

proc freq data=flag_pb_by_event noprint;
   tables flag_junction_annotated*flag_intron_retention*flag_pb_has_refseq*
          flag_event_redundant*flag_pb_in_refseq_list / out=event_pb_flags;
run;

proc print data=event_pb_flags;
run;


/*
                                   flag_pb_                    flag_pb_
 flag_junction_    flag_intron_      has_      flag_event_    in_refseq_
    annotated        retention      refseq      redundant        list       COUNT

Redundant events:
        0                1             0            1              .           20
        0                0             1            1              0            3
        0                1             1            1              0           45

Annotated junctions
        1                0             0            0              .         3555 ** <- limitation
        1                0             1            0              0           82 ** <- limitation
        1                0             1            0              1        57326 ** <- good

Unannotated junctions:
        0                0             0            0              .          344
        0                0             1            0              0           28
        0                0             1            0              1         1027

Intron retention:
        0                1             0            0              .          977
        0                1             1            0              0          313

*/

/* FRAGMENT */
