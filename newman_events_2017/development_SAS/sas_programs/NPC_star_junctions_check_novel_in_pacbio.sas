
/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* For "novel" junctions not in the RefSeq junction catalog, flag and count the number
   that are seen in the PacBio data */


data pb_junc;
   set evspl.splicing_events_annotations;
   where flag_intron_retention=0;
   intron_start=feature1_stop + 1;
   intron_stop=feature2_start;
   keep event_id gene_id chr transcript_id intron_start intron_stop flag_junction_annotated;
   rename event_id=pb_event_id transcript_id=pb_transcript_id gene_id=pb_gene_id
    flag_junction_annotated=flag_pb_junc_annot;
run;

data novel_junc;
  set event.NPC_star_junc_w_annot;
  where flag_in_catalog=0;
run;

proc sort data=pb_junc;
   by chr intron_start intron_stop;
proc sort data=novel_junc;
   by chr intron_start intron_stop;
run;

data novel_junc_pb_annot;
  merge novel_junc (in=in1) pb_junc (in=in2);
  by chr intron_start intron_stop;
  if in2 then flag_in_pacbio=1; else flag_in_pacbio=0;
  if in1 then output;
run;

/* Count: how many novel junctions are in the Pacbio data? */

proc freq data=novel_junc_pb_annot noprint;
  tables flag_in_pacbio*flag_in_catalog*flag_depth_NSC1_ge5*flag_depth_NSC2_ge5 / out=novel2pb_cnt;
run;

proc print data=novel2pb_cnt;
run;

/*

   flag_in_    flag_in_    flag_depth_    flag_depth_
    pacbio      catalog      NSC1_ge5       NSC2_ge5     COUNT

Not detected, not in PacBio:
       0           0            0              0         32059
Detected in one sample, not in PacBio:
       0           0            0              1            81
       0           0            1              0           387
Detected in both samples, not in PacBio:
       0           0            1              1           241

Not detected, in PacBio:
       1           0            0              0           801
Detected in one sample, in PacBio:
       1           0            0              1            17
       1           0            1              0            89
Detected in both samples, in PacBio:
       1           0            1              1            99
*/

data novel_junc_no_annot;
  set novel_junc_pb_annot;
  where flag_in_pacbio=0 and flag_depth_NSC1_ge5=1 and flag_depth_NSC2_ge5=1;
run;

/* Going to try remerging on either just intron_stop (as feature2_start) or intron_start (as feature1_start + 1) to see if perhaps
   the annotations are entirely correct... */

data novel_junc_no_annot_feat2;
  set novel_junc_no_annot;
  feature2_start=intron_stop;
run;

data novel_junc_no_annot_feat1;
  set novel_junc_no_annot;
  feature1_stop=intron_start-1;
run;


data pb_junc;
   set evspl.splicing_events_annotations;
   where flag_intron_retention=0;
   keep event_id gene_id chr transcript_id flag_junction_annotated  feature1_stop feature2_start;
   rename event_id=pb_event_id transcript_id=pb_transcript_id gene_id=pb_gene_id
    flag_junction_annotated=flag_pb_junc_annot;
run;

proc sort data=novel_junc_no_annot_feat2;
  by chr feature2_start;
proc sort data=pb_junc;
  by chr feature2_start;
run;

data novel_junc_feat2_check;
  merge novel_junc_no_annot_feat2 (in=in1) pb_junc (in=in2);
  by chr feature2_start;
  if in1;
run;


data novel_junc_feat2_check2; 
  set novel_junc_feat2_check;
  if intron_start - 1 = feature1_stop then flag_pb_match=1;
  else flag_pb_match=0;
run;

proc freq data=novel_junc_feat2_check2;
   tables flag_pb_match;
run;
           
/*

                             The FREQ Procedure

                                               Cumulative    Cumulative
     flag_pb_match    Frequency     Percent     Frequency      Percent
     ------------------------------------------------------------------
                 0         431      100.00           431       100.00

*/

 
proc sort data=novel_junc_no_annot_feat1;
  by chr feature1_stop;
proc sort data=pb_junc;
  by chr feature1_stop;
run;

data novel_junc_feat1_check;
  merge novel_junc_no_annot_feat1 (in=in1) pb_junc (in=in2);
  by chr feature1_stop;
  if in1;
run;


data novel_junc_feat1_check2; 
  set novel_junc_feat1_check;
  if intron_stop = feature2_start then flag_pb_match=1;
  else flag_pb_match=0;
run;

proc freq data=novel_junc_feat1_check2;
   tables flag_pb_match;
run;

/*


                                         Cumulative    Cumulative
lag_pb_match    Frequency     Percent     Frequency      Percent
-----------------------------------------------------------------
           0         462      100.00           462       100.00



Okay, so there are 241 junctions of high-confidence that are not obtainable
with either PacBio or via the event catalog

*/

