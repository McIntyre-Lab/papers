libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname event '!MCLAB/event_analysis/sas_data';


data pb_junc;
   set conesa.splicing_annotations_fixed;
   where transcript_id ne "";
   keep event_id event_size gene_id chr strand feature1_start feature1_stop feature2_start feature2_stop transcript_id;
   rename event_id=pb_junction_id;
run;

proc import datafile="!MCLAB/conesa_pacbio/created_files/pb_junc_all.tsv" out=pb_seq
    dbms=tab replace; guessingrows=685163;
run;

data pb_seq2;
  set pb_seq;
    sequence=tranwrd(sequence,"a","A");
  sequence=tranwrd(sequence,"c","C");
  sequence=tranwrd(sequence,"g","G");
  sequence=tranwrd(sequence,"t","T");
  sequence=tranwrd(sequence,"n","N");
run;


proc sort data=pb_seq2;
  by pb_junction_id;
proc sort data=pb_junc;
  by pb_junction_id;
run;

data pb_junc2seq;
   merge pb_junc (in=in1) pb_seq2 (in=in2);
   by pb_junction_id;
   if in1 and in2;
   rename event_size=pb_event_size
          gene_id=pb_gene_id;
run;

/* Get RefSeq junctions and sequences */

data rs_junc;
   set evspl.splicing_events_annot_refseq;
   if transcript_id="" then flag_junction_annotated=0; else flag_junction_annotated=1;
   keep event_id event_size gene_id chr strand feature1_start feature1_stop feature2_start feature2_stop flag_junction_annotated flag_intron_retention;
run;

data rs_seq;
  set evspl.flag_event_redundant_seq;
  rename sequence=event_sequence;
run;

proc sort data=rs_junc;
  by event_id;
proc sort data=rs_seq;
  by event_id;
run;

data rs_junc2seq;
    merge rs_junc (in=in1) rs_seq (in=in2);
   by event_id;
   if in1 and in2;
run;

data gene2keep;
  set event.feature2xs2gene_nomulti;
  keep gene_id;
run;

proc sort data=rs_junc2seq;
  by gene_id;
proc sort data=gene2keep nodup;
  by gene_id;
run;

data rs_junc2seq_nomult;
  merge rs_junc2seq (in=in1) gene2keep (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=rs_junc2seq_nomult;
   by chr feature1_start feature1_stop feature2_start feature2_stop strand;
proc sort data=pb_junc2seq;
   by chr feature1_start feature1_stop feature2_start feature2_stop strand;
run;

data pb2rs_junc2seq_nomulti nopb;
   merge rs_junc2seq_nomult (in=in1) pb_junc2seq (in=in2);
   by chr feature1_start feature1_stop feature2_start feature2_stop strand;
   if in1 and in2 then do;
          flag_has_pb_match=1;
          output pb2rs_junc2seq_nomulti;
          end;
   else if in1 then do;
          flag_has_pb_match=0;
          output pb2rs_junc2seq_nomulti;
          end;
   else output nopb;
run;

proc freq data=pb2rs_junc2seq_nomulti;
   table flag_has_pb_match*flag_junction_annotated;
run;

/* 

  flag_has_pb_match
            flag_junction_annotated

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |1788449 | 148111 |1936560
           |  90.11 |   7.46 |  97.58
           |  92.35 |   7.65 |
           |  99.93 |  75.99 |
  ---------+--------+--------+
         1 |   1303 |  46803 |  48106
           |   0.07 |   2.36 |   2.42
           |   2.71 |  97.29 |
           |   0.07 |  24.01 |
  ---------+--------+--------+
  Total     1789752   194914  1984666
              90.18     9.82   100.00

Okay. I expect ~1K unannotated junctions to therefore match to PB transcripts!

Now, check if events have the same sequence (they should!)

*/

data check_seq;
   set pb2rs_junc2seq_nomulti;
   if flag_has_pb_match=1 then do;
   if event_sequence=sequence then flag_seq_ok=1;
   else flag_seq_ok=0;
   if event_size=pb_event_size then flag_size_ok=1;
   else flag_size_ok=0; *this should be the same as flag_seq_ok;
   end;
run;

proc freq data=check_seq;
   tables flag_seq_ok*flag_size_ok;
run;

proc freq data=check_seq;
   tables flag_seq_ok*flag_junction_annotated;
run;

/*
     flag_seq_ok
               flag_junction_annotated

     Frequency|
     Percent  |
     Row Pct  |
     Col Pct  |       0|       1|  Total
     ---------+--------+--------+
            0 |     93 |   1168 |   1261
              |   0.19 |   2.43 |   2.62
              |   7.38 |  92.62 |
              |   7.14 |   2.50 |
     ---------+--------+--------+
            1 |   1210 |  45635 |  46845
              |   2.52 |  94.86 |  97.38
              |   2.58 |  97.42 |
              |  92.86 |  97.50 |
     ---------+--------+--------+
     Total        1303    46803    48106
                  2.71    97.29   100.00

So sequences generally match. Now to check which are actually detected
*/

data event_on;
  set event.splicing_on_apn_gt0;
  keep event_id flag_splicing_on;
run;

proc sort data=event_on;
  by event_id;
proc sort data=check_seq;
  by event_id;
run;

data check_on;
  merge check_seq (in=in1) event_on (in=in2);
  by event_id;
  if in1 and in2;
run;

proc freq data=check_on noprint;
  tables flag_junction_annotated*flag_has_pb_match*flag_seq_ok*flag_size_ok*flag_splicing_on*flag_redundant_sequence / out=all_pb2rs_checks;
run;

proc print data=all_pb2rs_checks;
run;


/*

                  flag_                         flag_
flag_junction_   has_pb_    flag_    flag_    splicing_   flag_redundant_
   annotated      match    seq_ok   size_ok       on          sequence        COUNT

       0            0         .        .          0              0          1734680
       0            0         .        .          0              1            16401
       0            0         .        .          1              0            37368
       0            1         0        1          0              0               18
       0            1         0        1          0              1                3
       0            1         0        1          1              0               72
       0            1         1        1          0              0              193
       0            1         1        1          0              1                4
       0            1         1        1          1              0             1013
       1            0         .        .          0              0           102778
       1            0         .        .          0              1             2322
       1            0         .        .          1              0            43011
       1            1         0        1          0              0               61
       1            1         0        1          0              1               28
       1            1         0        1          1              0             1079
       1            1         1        1          0              0             1629
       1            1         1        1          0              1               91
       1            1         1        1          1              0            43915


Again, there SHOULD BE ~1000 unannotated junctions that match to PB !

*/

/* Make permenant -- can use to check later */

data event.pacbio2refseq_junction_check;
   set check_on;
   keep event_id pb_junction_id transcript_id pb_gene_id gene_id flag_junction_annotated flag_intron_retention flag_has_pb_match flag_seq_ok flag_size_ok flag_splicing_on flag_redundant_sequence
        event_sequence sequence;
run;








