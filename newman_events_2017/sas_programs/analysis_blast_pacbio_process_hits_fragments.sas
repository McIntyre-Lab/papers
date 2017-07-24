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


/* FRAGMENT */

*non-redundant fragments;
data nr_frags;
  set event.flag_redundant_frag_rfsq_blast;
  if flag_blast_other=1 then flag_frag_redundant=1;
  else flag_frag_redundant=0;
  keep fragment_id flag_frag_redundant flag_frag_on flag_from_exp_gene;
run;

*BLAST hits;

data annot_hits;
   set event.frag_pacbio_blast_hits_known;
   flag_annot_hit=1;
   flag_unannot_hit=0;
   flag_novel_hit=0;
run;

data unannot_hits;
   set event.frag_pacbio_blast_hits_other;
   flag_annot_hit=0;
   flag_unannot_hit=1;
   flag_novel_hit=0;
run;

data novel_hits;
   set event.frag_pacbio_blast_hits_no_refseq;
   flag_annot_hit=0;
   flag_unannot_hit=0;
   flag_novel_hit=1;
run;

data blast_hits;
  set annot_hits unannot_hits novel_hits;
run;

proc sort data=blast_hits;
  by fragment_id;
proc sort data=nr_frags nodup;
  by fragment_id;
run;

data blast_hits_pb2rs_nr;
  merge blast_hits (in=in1) nr_frags (in=in2);
  by fragment_id;
  if in1 and in2 then output;
  else if in1 then do;
     flag_frag_redundant=0;
     output;
     end;
run;

*check;

proc freq data=blast_hits_pb2rs_nr;
  tables flag_has_refseq_id*flag_frag_redundant;
run;

/*
 flag_has_refseq_id
           flag_frag_redundant

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  28901 |    717 |  29618
          |  36.75 |   0.91 |  37.66
          |  97.58 |   2.42 |
          |  37.40 |  52.37 |
 ---------+--------+--------+
        1 |  48370 |    652 |  49022
          |  61.51 |   0.83 |  62.34
          |  98.67 |   1.33 |
          |  62.60 |  47.63 |
 ---------+--------+--------+
 Total       77271     1369    78640
             98.26     1.74   100.00


*/

/* Flag if PB transcript is in RefSeq list */

data flag_pb_in_refseq;
   set blast_hits_pb2rs_nr;
   if flag_has_refseq_id=1 then do;
   if count(transcript_id_cat,strip(transcript_id)) > 0 then flag_pb_in_refseq_list=1;
   else flag_pb_in_refseq_list=0;
   end;
run;

/* See if event aligns to any RefSeq transcript */

proc sort data=flag_pb_in_refseq;
   by fragment_id;
run;

proc means data=flag_pb_in_refseq noprint;
   by fragment_id;
   var flag_has_refseq_id flag_annot_hit flag_unannot_hit flag_novel_hit flag_frag_redundant
        flag_pb_in_refseq_list  flag_frag_on flag_from_exp_gene;
   output out=flag_pb_by_frag max=;
run;

proc freq data=flag_pb_by_frag noprint;
   tables flag_has_refseq_id*flag_annot_hit*flag_unannot_hit*flag_novel_hit*
          flag_frag_redundant*flag_pb_in_refseq_list*flag_frag_on*flag_from_exp_gene / out=frag_pb_flags;
run;

proc print data=frag_pb_flags;
run;


/*

            flag_    flag_    flag_               flag_pb_              flag_
flag_has_  annot_  unannot_  novel_  flag_frag_  in_refseq_   flag_     from_
refseq_id    hit      hit      hit    redundant     list     frag_on  exp_gene  COUNT  PERCENT

These I don't really care about:
    0         0        0        1         0           .         .         .       299    .
    1         0        1        0         0           0         .         .        91    .
    1         1        0        0         0           1         .         .       990    .
    1         1        0        1         0           1         .         .       561    .
    1         1        1        0         0           1         .         .        14    .
    1         1        1        1         0           1         .         .        5    .
    1         0        1        1         0           0         .         .        86    .

These I don't really care about:
    0         0        0        1         1           .         0         1         5    .
    1         0        1        0         1           0         0         1         1   0.0034
    1         0        1        1         1           0         0         1         2   0.0068
    1         1        0        1         1           1         0         1         1   0.0034
    0         0        0        1         0           .         0         1        32    .
    1         0        1        0         0           0         0         1         2   0.0068
    1         1        0        0         0           1         0         1       220   0.7488
    1         1        0        1         0           1         0         1        16   0.0545

These I don't really care about (redundant fragments, but can report:
    1         1        1        0         1           1         1         1        5   0.01702
    0         0        0        1         1           .         1         1        43    .
    1         0        1        0         1           0         1         1        56   0.1906
    1         0        1        1         1           0         1         1        56   0.1906
    1         1        0        0         1           1         1         1       102   0.3472
    1         1        0        1         1           1         1         1        95   0.3233
    1         1        1        1         1           1         1         1       56   0.19061

These ones I care about:
    0         0        0        1         0           .         1         1      1857    .
    1         0        1        0         0           0         1         1        48   0.1634
    1         0        1        1         0           0         1         1        15   0.0511
    1         1        0        0         0           1         1         1     17897  60.9156
    1         1        0        1         0           1         1         1     10748  36.5827
    1         1        1        0         0           1         1         1       29   0.09871
    1         1        1        1         0           1         1         1       31   0.10551

1857+48+15+47897+10748+29+31 expressed, non-redundant fragments in expressed genes have hits
*/


data event.blast_pacbio_fragments_in_refseq;
  set flag_pb_in_refseq;
run;

