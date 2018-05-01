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
proc sort data=event.simul_star_junctions;
  by chr intron_start intron_stop;
run;

data star_junc_w_annot;
  merge event.simul_star_junctions (in=in1) junc2intron (in=in2);
  by chr intron_start intron_stop;
  if in2 then flag_in_catalog=1; else flag_in_catalog=0;
  if in1 then output;
run;

/* Check : how many junctions are "annotated" in STAR output, vs number annotated in catalog ? */

proc freq data=star_junc_w_annot;
  tables flag_junction_annotated*flag_junc_annotated_cat;
run;

/*
 flag_junction_annotated
           flag_junc_annotated_cat

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 | 188103 |  56213 | 244316
          |  49.93 |  14.92 |  64.85
          |  76.99 |  23.01 |
          |  96.22 |  31.01 |
 ---------+--------+--------+
        1 |   7398 | 125036 | 132434
          |   1.96 |  33.19 |  35.15
          |   5.59 |  94.41 |
          |   3.78 |  68.99 |
 ---------+--------+--------+
 Total      195501   181249   376750
             51.89    48.11   100.00
Frequency Missing = 109080

Of the total number (no filtering):

109080 junctions not in the catalog
188103 junctions called as unannotated by both STAR and the junction catalog
56213 junctions called as annotated in the junction catalog, but not by STAR
7398 junctions called as annotated by STAR, but are unannotated in the junction catalog
125036 that are called as annotated in both

*/

/* Look at one that are discordant */

data junc_check;
  set star_junc_w_annot;
  where (flag_junction_annotated=0 and flag_junc_annotated_cat=1)
  or (flag_junction_annotated=1 and flag_junc_annotated_cat=0);
run;

/* 
Column might only be useful if you supply a list of annotated junctions...
*/

/* Check for each test: how many junctions with at least 5 reads are not in the junction catalog? */

proc freq data=star_junc_w_annot noprint;
   tables flag_depth_sim1_test1_ge5*flag_depth_sim2_test1_ge5*
          flag_depth_sim3_test1_ge5*flag_in_catalog*flag_junc_annotated_cat / out=test1_junc_count;
run;

proc freq data=star_junc_w_annot noprint;
   tables flag_depth_sim1_test2_ge5*flag_depth_sim2_test2_ge5*
          flag_depth_sim3_test2_ge5*flag_in_catalog*flag_junc_annotated_cat / out=test2_junc_count;
run;

proc print data=test1_junc_count;
run;

/*
  sim1_test1_    sim2_test1_    sim3_test1_    flag_in_    annotated_
      ge5            ge5            ge5         catalog        cat        COUNT

Not detected:
       0              0              0             0            .         83137
       0              0              0             1            0        156521
       0              0              0             1            1         58961
One test data only:
       0              0              1             0            .          7666
       0              0              1             1            0         10728
       0              0              1             1            1         13922

       0              1              0             0            .          7741
       0              1              0             1            0         11299
       0              1              0             1            1         15180

       1              0              0             0            .          7683
       1              0              0             1            0         11135
       1              0              0             1            1         14499

Two test data:
       0              1              1             0            .           891
       0              1              1             1            0          1675
       0              1              1             1            1         15782
       1              0              1             0            .           802
       1              0              1             1            0          1614
       1              0              1             1            1         14503
       1              1              0             0            .           935
       1              1              0             1            0          1645
       1              1              0             1            1         14538

All test data:
       1              1              1             0            .           225
       1              1              1             1            0           884
       1              1              1             1            1         33864

Looks like when we only pick out high-confidence junctions (ie, we see them in most of our sample)
then the ones that are detected are almost all annotated junctions

Calculate proportions for each

*/

proc print data=test2_junc_count;
run;

/*

 flag_depth_    flag_depth_    flag_depth_                flag_junc_
 sim1_test2_    sim2_test2_    sim3_test2_    flag_in_    annotated_
     ge5            ge5            ge5         catalog        cat        COUNT

Not detected:
      0              0              0             0            .         83824
      0              0              0             1            0        141610
      0              0              0             1            1         63637

One test data only:
      0              0              1             0            .          8129
      0              0              1             1            0         16232
      0              0              1             1            1         15041
      0              1              0             0            .          7385
      0              1              0             1            0         14873
      0              1              0             1            1         13368
      1              0              0             0            .          7315
      1              0              0             1            0         14926
      1              0              0             1            1         14390

Two test data:
      0              1              1             0            .           762
      0              1              1             1            0          2382
      0              1              1             1            1         14744
      1              0              1             0            .           713
      1              0              1             1            0          2248
      1              0              1             1            1         14418
      1              1              0             0            .           763
      1              1              0             1            0          2296
      1              1              0             1            1         14432

All test data:
      1              1              1             0            .           189
      1              1              1             1            0           934
      1              1              1             1            1         31219

Same story as above.
*/

data event.simul_star_junc_w_annot;
  set  star_junc_w_annot;
run;

