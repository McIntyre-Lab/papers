/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Count STAR junctions also in junction catalog */

data junctions;
   set evspl.splicing_events_annotations;
   keep event_id gene_id chr strand transcript_id flag_junction_annotated feature1_stop feature2_start;
run;

data junc2intron;
  set junctions;
  where event_id ^? "intron";
  intron_start=feature1_stop + 1;
  intron_stop=feature2_start;
  rename strand=event_strand flag_junction_annotated=flag_junc_annotated_cat;
run;

proc sort data=junc2intron;
  by chr intron_start intron_stop;
proc sort data=event.NPC_star_junctions;
  by chr intron_start intron_stop;
run;

data star_junc_w_annot;
  merge event.NPC_star_junctions (in=in1) junc2intron (in=in2);
  by chr intron_start intron_stop;
  if in2 then flag_in_catalog=1; else flag_in_catalog=0;
  if in1 then output;
run;

/* Check for each test: how many junctions with at least 5 reads are not in the junction catalog? */

proc freq data=star_junc_w_annot noprint;
   tables flag_depth_NSC1_ge5*flag_depth_NSC2_ge5*flag_in_catalog*flag_junc_annotated_cat / out=NPC_junc_count;
run;

proc print data=NPC_junc_count;
run;

/*

                                           flag_junc_
 flag_depth_    flag_depth_    flag_in_    annotated_
   NSC1_ge5       NSC2_ge5      catalog        cat       COUNT
Not detected:
      0              0             0            .        75151
      0              0             1            0         4946
      0              0             1            1        14151
Single sample only:
      0              1             0            .         1492
      0              1             1            0           94
      0              1             1            1         1567
      1              0             0            .         5865
      1              0             1            0          346
      1              0             1            1         7189
Both samples
      1              1             0            .         8044
      1              1             1            0          977
      1              1             1            1        42358

Most detected junctions are in the PacBio data, however there is a substantial number of junctions
that are NOT seen in the PacBio catalog

Suspect most of those "novel" are from known transcripts not in the PacBio data...
*/

data event.NPC_star_junc_w_annot_pacbio;
  set  star_junc_w_annot;
run;

