ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* REFSEQ BLAST PROCESSING:

 1. If event aligned to non-originating transcripts, then flag as ambiguous and remove from further analysis
2. Else event is okay if they align to at least one of their originating transcripts.
3. Fragments are okay if they align to at least one of their originating transcripts.
	- still flag if it goes elsewhere but don't remove yet!
*/

/* Events */

data event_redundant;
   set event.refseq_blast_redundancy_class;
   keep event_id;
run;

data event_blast;
   set event.refseq_blast_hits_w_annot;
run;

proc sort data=event_redundant nodup;
  by event_id; *363 redundant events -- not too many!;
proc sort data=event_blast nodup;
  by event_id;
run; *127266 events total;

data event_blast_flag_redundant oops;
   merge event_blast (in=in1) event_redundant (in=in2);
   by event_id;
   if in1 and in2 then do;
         flag_event_redundant=1;
         output event_blast_flag_redundant;
         end;
   else if in1 then do;
         flag_event_redundant=0;
         output event_blast_flag_redundant;
         end;
   else output oops;
run;

proc freq data=event_blast_flag_redundant noprint;
   tables flag_junction_annotated*flag_intron_retention*flag_annotated_blast_hit*
          flag_unannotated_blast_hit*flag_event_redundant / out=event_redund_count;
run;

proc print data=event_redund_count;
run;

/*

flag_junction_  flag_intron_  flag_annotated_  flag_unannotated_  flag_event_
   annotated      retention      blast_hit         blast_hit       redundant    COUNT  PERCENT

       0              0              0                 1               1           20   0.0157
       0              1              0                 1               1          343   0.2695
       1              0              0                 1               0            1   0.0008
       1              0              1                 0               0       126137  99.1129
       1              0              1                 1               0          765   0.6011

All IR events with hits to RefSeq are redundant sequence
20 unannotated junctions that have redundant sequence
The rest are annotated junctions

The 1 that is annotated but all hits are "unannotated" is apparently due to a bug in my code
somewhere. It in fact matches to all its actual transcripts

*/

/* Fragments */

data frag_blast_known;
   set event.frag_refseq_blast_hits_known;
   keep fragment_id;
run;

data frag_blast_other;
   set event.frag_refseq_blast_hits_other;
   keep fragment_id;
run;

proc sort data=frag_blast_known nodup;
   by fragment_id;
proc sort data=frag_blast_other nodup;
   by fragment_id;
run;

data frag_blast_flags;
   merge frag_blast_known (in=in1) frag_blast_other (in=in2);
   by fragment_id;
   if in1 then flag_blast_known=1; else flag_blast_known=0;
   if in2 then flag_blast_other=1; else flag_blast_other=0;
run;

proc freq data=frag_blast_flags;
   tables flag_blast_known*flag_blast_other;
run;

/*
  flag_blast_known
            flag_blast_other

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |    299 |    299
           |   0.00 |   0.19 |   0.19
           |   0.00 | 100.00 |
           |   0.00 |  12.13 |
  ---------+--------+--------+
         1 | 151718 |   2166 | 153884
           |  98.40 |   1.40 |  99.81
           |  98.59 |   1.41 |
           | 100.00 |  87.87 |
  ---------+--------+--------+
  Total      151718     2465   154183
              98.40     1.60   100.00
*/

/* Flag if event is in expressed gene, and is detected then count only these ones */

data frag_exp_gene;
   set event.feature2xs2gene;
   keep feature_id;
   rename feature_id=fragment_id;
run;

data frag_on;
   set event.fragments_on_apn_gt0;
   where flag_fragment_on=1;
   keep fragment_id;
run;

proc sort data=frag_exp_gene nodup;
   by fragment_id;
proc sort data=frag_on nodup;
   by fragment_id;
run;

data frag_blast_flags_exp_on;
  merge frag_blast_flags (in=in1) frag_exp_gene (in=in2) frag_on (in=in3);
  by fragment_id;
  if in2 then flag_from_exp_gene=1; else flag_from_exp_gene=0;
  if in3 then flag_frag_on=1; else flag_frag_on=0;
  if in1 then output;
run;


proc freq data=frag_blast_flags_exp_on;
   where flag_from_exp_gene=1;
   tables flag_blast_known*flag_blast_other;
run;

/*
 flag_blast_known
           flag_blast_other

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |      0 |    298 |    298
          |   0.00 |   0.19 |   0.19
          |   0.00 | 100.00 |
          |   0.00 |  12.09 |
 ---------+--------+--------+
        1 | 151718 |   2166 | 153884
          |  98.40 |   1.40 |  99.81
          |  98.59 |   1.41 |
          | 100.00 |  87.91 |
 ---------+--------+--------+
 Total      151718     2464   154182
             98.40     1.60   100.00

*/

proc freq data=frag_blast_flags_exp_on;
   where flag_frag_on=1;
   tables flag_blast_known*flag_blast_other;
run;

/*
  flag_blast_known
            flag_blast_other

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |     90 |     90
           |   0.00 |   0.09 |   0.09
           |   0.00 | 100.00 |
           |   0.00 |   9.68 |
  ---------+--------+--------+
         1 |  99795 |    840 | 100635
           |  99.08 |   0.83 |  99.91
           |  99.17 |   0.83 |
           | 100.00 |  90.32 |
  ---------+--------+--------+
  Total       99795      930   100725
              99.08     0.92   100.00

*/

proc freq data=frag_blast_flags_exp_on;
   where flag_from_exp_gene=1 and flag_frag_on=1;
   tables flag_blast_known*flag_blast_other;
run;

/*
  flag_blast_known
            flag_blast_other

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |      0 |     89 |     89
           |   0.00 |   0.09 |   0.09
           |   0.00 | 100.00 |
           |   0.00 |   9.58 |
  ---------+--------+--------+
         1 |  99795 |    840 | 100635
           |  99.08 |   0.83 |  99.91
           |  99.17 |   0.83 |
           | 100.00 |  90.42 |
  ---------+--------+--------+
  Total       99795      929   100724
              99.08     0.92   100.00
*/


/* Decision: if ambig with no "known" hit, then flag as redundant and drop */

data event.flag_redundant_event_rfsq_blast;
   set event_blast_flag_redundant;
run;

data event.flag_redundant_frag_rfsq_blast;
   set frag_blast_flags_exp_on;
run;

