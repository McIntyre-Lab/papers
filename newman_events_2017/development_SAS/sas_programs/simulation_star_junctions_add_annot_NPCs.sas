/* Libraries */

ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Count STAR junctions also in junction catalog */

data junctions;
   set evspl.splicing_events_annot_refseq;
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
     0              0             0            .        32804
     0              0             1            0         7218
     0              0             1            1        54954
Single sample only:

     0              1             0            .           98
     0              1             1            0          107
     0              1             1            1         2956
     1              0             0            .          473
     1              0             1            0          524
     1              0             1            1        12401
Both samples
     1              1             0            .          332
     1              1             1            0         1137
     1              1             1            1        49413

So similar to the simulated data, most detectable junctions are in the junction catalog

Most of the "not detected" (i.e. has very low coverage) junctions are novel junctions
*/

data event.NPC_star_junc_w_annot;
  set  star_junc_w_annot;
run;

