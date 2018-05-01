/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the 10000 transcript simulation, I want a table that compared Event Analysis to iReckon,
   with the following info:
   (a) Number of correctly-identified transcripts
   (b) Number of related Refseq transcripts (i.e. from the same gene as simulated transcripts)
       that are detected
   (c) Number of related non-Refseq transcripts (iReckon only)
   (d) Number of unrelated transcripts (i.e. from other genes)
   (e) Number of missing transcripts (simulated but not correctly identified

For Events, I will use to criteria:
100% events detected at APN>0
>=75% events detected at APN>=5

*/

/* (1) Get simulated transcripts and their related RefSeq transcripts */

data xs10000;
  set event.polyester_xs_list_10k;
run;

data xs2gene;
   set event.feature2xs2gene;
   drop feature_id;
run;

proc sort data=xs10000;
  by transcript_id;
proc sort data=xs2gene nodup;
  by transcript_id;
run;

data rs_gene;
  merge xs10000 (in=in1) xs2gene (in=in2);
  by transcript_id;
  if in1 and in2;
  keep gene_id;
run;

proc sort data=rs_gene nodup;
  by gene_id; *7678 genes;
proc sort data=xs2gene;
  by gene_id;
run;

data rs_xscript;
  merge rs_gene (in=in1) xs2gene (in=in2);
  by gene_id;
  if in1 and in2;
  keep transcript_id;
run; *43198 RefSeq transcripts;

proc sort data=xs10000;
   by transcript_id;
proc sort data=rs_xscript;
   by transcript_id;
run;

data related_refseq;
  merge rs_xscript (in=in1) xs10000 (in=in2);
  by transcript_id;
  if in2 then delete;
run; *33198 possible related RefSeq transcripts;

/* Events (RR1, RR2) - what transcripts are detected, what are missing? */

data event100;
   set event.bin_xs_by_dtct_apn0_10k;
   where perc_features_dtct = 1;
   keep transcript_id;
run;

data event75;
   set event.bin_xs_by_dtct_apn0_10k;
   where perc_features_dtct >= 0.75;
   keep transcript_id;
run;


data event75_apn5;
   set event.bin_xs_by_dtct_apn5_10k;
   where perc_features_dtct >= 0.75 ;
   keep transcript_id;
run;

proc sort data=event100;
  by transcript_id;
proc sort data=event75;
  by transcript_id;
proc sort data=event75_apn5;
  by transcript_id;
proc sort data=related_refseq;
  by transcript_id;
proc sort data=xs10000;
  by transcript_id;
run;

/* Events part of table */

data event_summary_100;
  merge xs10000 (in=in1) related_refseq (in=in2) event100 (in=in3) ;
  by transcript_id;
  if in1 then flag_sim10000=1; else flag_sim10000=0;
  if in2 then flag_related_refseq=1; else flag_related_refseq=0;
  if in3 then flag_in_event100_apn0=1; else flag_in_event100_apn0=0;
run;

data event_summary_75;
  merge xs10000 (in=in1) related_refseq (in=in2) event75 (in=in3);
  by transcript_id;
  if in1 then flag_sim10000=1; else flag_sim10000=0;
  if in2 then flag_related_refseq=1; else flag_related_refseq=0;
  if in3 then flag_in_event75_apn0=1; else flag_in_event75_apn0=0;
run;



data event_summary_75_5;
  merge xs10000 (in=in1) related_refseq (in=in2) event75_apn5 (in=in3);
  by transcript_id;
  if in1 then flag_sim10000=1; else flag_sim10000=0;
  if in2 then flag_related_refseq=1; else flag_related_refseq=0;
  if in3 then flag_in_event75_apn5=1; else flag_in_event75_apn5=0;
run;

proc freq data=event_summary_100 noprint;
  tables flag_sim10000*flag_related_refseq*flag_in_event100_apn0 / out=ea_100;
run;

proc freq data=event_summary_75 noprint;
  tables flag_sim10000*flag_related_refseq*flag_in_event75_apn0 / out=ea_75;
run;


proc freq data=event_summary_75_5 noprint;
  tables flag_sim10000*flag_related_refseq*flag_in_event75_apn5 / out=ea_75_5;
run;

proc print data=ea_100;
run;
proc print data=ea_75;
run;
proc print data=ea_75_5;
run;

/* EVENTS: 100%, APN>0
                                flag_in_
    flag_     flag_related_    event100_
  sim10000        refseq          apn0      COUNT

      0             0              1          343
      0             1              0        28083
      0             1              1         5115
      1             0              0          728
      1             0              1         9272

# Correct			9272
# Related RefSeq	5115
# Related other		0
# Unrelated			343
# Missing			728

EVENTS: 75+%, APN>0

                             flag_in_
  flag_     flag_related_    event75_
sim10000        refseq         apn0      COUNT

    0             0              1         399
    0             1              0        8229
    0             1              1       24969
    1             0              0          63
    1             0              1        9937

# Correct			9937
# Related RefSeq	24969
# Related other		0
# Unrelated			399
# Missing			63


EVENTS: 75+%, APN>=5
                              flag_in_
   flag_     flag_related_    event75_
 sim10000        refseq         apn5      COUNT

     0             0              1         134
     0             1              0        9131
     0             1              1       24067
     1             0              0         119
     1             0              1        9881

# Correct			9881
# Related RefSeq	24067
# Related other		0
# Unrelated			134
# Missing			119

*/


/* iReckon: I am going to use the "union" iReckon transcriptome here.
   I am not going to use the iReckon-to-Refseq BLAST results here, just those IDs
   that iReckon correctly assigns

*/

data ir;
  set event.polyester_10k_simulation_ireckon;
run;

/*   Because of the way that iReckon assigns gene IDs when none are provided (yet, I can't find
   info on this in the IR documentation...), I am going to assign RefSeq gene IDs to the iReckon IDs.
   Then assign the Refseq gene ID to each iR transcript with the same iR gene ID */

data ir_genes;
  set ir;
  length refseq_id $100.;
  if index(gene_id,"NM_") > 0 or index(gene_id,"NR_") > 0
     or index(gene_id,"XM_") > 0 or index(gene_id,"XR_") > 0
     then refseq_id=tranwrd(gene_id,"IntronRetention","");
  else refseq_id=tranwrd(transcript_id,"IntronRetention","");
  keep sample_id gene_id transcript_id refseq_id strand;
run;

data xs2gene;
   set event.feature2xs2gene;
   keep gene_id transcript_id;
   rename transcript_id=refseq_id gene_id=refseq_gene_id;
run;

proc sort data=ir_genes;
  by refseq_id;
proc sort data=xs2gene nodup;
  by refseq_id;
run;


data ir_gene_w_refseq;
  merge ir_genes (in=in1) xs2gene (in=in2);
  by refseq_id;
  if in1 and in2;
  keep sample_id gene_id refseq_gene_id  strand;
run;

proc sort data=ir_gene_w_refseq nodup;
  by sample_id  strand gene_id refseq_gene_id;
run;

proc sort data=ir;
  by sample_id  strand gene_id;
run;

data ir_w_refseq_gene_id;
  merge ir (in=in1) ir_gene_w_refseq (in=in2);
  by sample_id strand gene_id;
  if in2 then flag_refseq_gene=1; else flag_refseq_gene=0;
  if in1 then output;
run;

/* Collapse transcripts -- some transcripts have the same RefSeq ID, but different gene_ids
   so I am going to */

data ir_w_refseq_gene_id2;
  set ir_w_refseq_gene_id;
  length transcript_id2 $100.;
  if index(transcript_id,"NM_") > 0 or index(transcript_id,"NR_") > 0 or 
     index(transcript_id,"XM_") > 0 or index(transcript_id,"XR_") > 0 
     then transcript_id2=transcript_id;
  else transcript_id2=catx("|",gene_id, transcript_id);
  run;

proc sort data=ir_w_refseq_gene_id2;
  by transcript_id2 ;
run;

proc means data=ir_w_refseq_gene_id2 noprint;
  by transcript_id2;
  var flag_refseq_gene;
  output out=ir_collapse max=;
run;

data xs10000;
  set event.polyester_xs_list_10k;
  rename transcript_id=transcript_id2;
run;

proc sort data=xs10000;
  by transcript_id2;
proc sort data=ir_collapse;
  by transcript_id2;
run;

data ir_collapse_flag_10k;
  merge ir_collapse (in=in1) xs10000 (in=in2);
  by transcript_id2;
  if in1 then flag_detected=1; else flag_detected=0;
  if in2 then flag_simulated=1; else flag_simulated=0;
  if (index(transcript_id2,"NM_") > 0 or index(transcript_id2,"NR_") > 0 or 
     index(transcript_id2,"XM_") > 0 or index(transcript_id2,"XR_") > 0 )
     and index(transcript_id2,"Intron") = 0 and index(transcript_id2,"novel") = 0
          and index(transcript_id2,"unspliced") = 0
     then flag_refseq_xs=1; else flag_refseq_xs=0;
run;

proc freq data=ir_collapse_flag_10k noprint;
  tables flag_detected*flag_simulated*flag_refseq_gene*flag_refseq_xs / out=ireckon_10k_count;
run;

proc print data=ireckon_10k_count;
run;

/*
                           flag_      flag_
   flag_       flag_      refseq_    refseq_
 detected    simulated      gene        xs      COUNT    PERCENT

     0           1           .          1        1942      .
     1           0           0          0        1127     6.3914
     1           0           1          0        5857    33.2161
     1           0           1          1        2591    14.6940
     1           1           1          1        8058    45.6984

CORRECT			8058	transcripts that were selected for simulation
RELATED_RS		4324	RefSeq transcripts that were also estimated, but not from simulation
RELATED_OTHER	4124	Non-RefSeq transcripts estimated and belong to same genes as a Refseq gene
UNRELATED		1127	All other transcripts
MISSING			1942	transcripts simulated that were not estimated
*/

data ir_collapse2;
  set ir_collapse;
  where _FREQ_ = 6;
run;


proc sort data=xs10000;
  by transcript_id2;
proc sort data=ir_collapse2;
  by transcript_id2;
run;

data ir_collapse_flag_10k;
  merge ir_collapse2 (in=in1) xs10000 (in=in2);
  by transcript_id2;
  if in1 then flag_detected=1; else flag_detected=0;
  if in2 then flag_simulated=1; else flag_simulated=0;
  if (index(transcript_id2,"NM_") > 0 or index(transcript_id2,"NR_") > 0 or 
     index(transcript_id2,"XM_") > 0 or index(transcript_id2,"XR_") > 0 )
     and index(transcript_id2,"Intron") = 0 and index(transcript_id2,"novel") = 0
          and index(transcript_id2,"unspliced") = 0
     then flag_refseq_xs=1; else flag_refseq_xs=0;
run;

proc freq data=ir_collapse_flag_10k noprint;
  tables flag_detected*flag_simulated*flag_refseq_gene*flag_refseq_xs / out=ireckon_10k_count;
run;

proc print data=ireckon_10k_count;
run;

/*

                           flag_      flag_
   flag_       flag_      refseq_    refseq_
 detected    simulated      gene        xs      COUNT    PERCENT

     0           1           .          1        3363      .
     1           0           0          0         228     2.6475
     1           0           1          0         502     5.8291
     1           0           1          1        1245    14.4566
     1           1           1          1        6637    77.0669



CORRECT			6637	transcripts that were selected for simulation
RELATED_RS		1245	RefSeq transcripts that were also estimated, but not from simulation
RELATED_OTHER	502		Non-RefSeq transcripts estimated and belong to same genes as a Refseq gene
UNRELATED		228		All other transcripts
MISSING			3363	transcripts simulated that were not estimated
*/

