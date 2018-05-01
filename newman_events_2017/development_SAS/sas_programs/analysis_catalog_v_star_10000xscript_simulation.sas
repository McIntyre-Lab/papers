/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the 10000 transcript simulation, import the STAR junction results and EA junction coverage counts and check
   to see which method captures the junctions of these genes */

/* Flag junctions from transcripts selected for simulation
   Flag junctions from gene with transcripts simulated
   Keep junction type (annot, NIC, border) */

data xscript_sim;
   set event.polyester_xs_list_10k;
run;

data event2xs;
   set event.feature2xs2gene;
   keep feature_id transcript_id;
   rename feature_id=event_id;
run;

proc sort data=xscript_sim;
   by transcript_id;
proc sort data=event2xs nodup;
   by transcript_id event_id;
run;

data event_sim;
  merge event2xs (in=in1) xscript_sim (in=in2);
  by transcript_id;
  if in1 and in2;
  keep event_id;
run;

data xs2gene;
   set event.feature2xs2gene;
   keep gene_id transcript_id;
run;

proc sort data=xscript_sim;
   by transcript_id;
proc sort data=xs2gene nodup;
   by transcript_id gene_id;
run;

data gene_w_sim_xs;
   merge xs2gene (in=in1) xscript_sim (in=in2);
   by transcript_id;
   if in1 and in2;
   keep gene_id;
run;

data event2gene;
   set evspl.splicing_events_annot_refseq;
   keep event_id gene_id;
run;

proc sort data=gene_w_sim_xs nodup;
  by gene_id;
proc sort data=event2gene;
  by gene_id;
run;

data event_gene_sim;
  merge event2gene (in=in1) gene_w_sim_xs (in=in2);
  by gene_id;
  if in1 and in2;
  keep event_id;
run;


data event2uniq;
  set eventloc.unique_junction2event_mm10;
  keep junction_id seq_name event_id;
run;

proc sort data=event2uniq;
   by event_id;
proc sort data=event_gene_sim nodup;
   by event_id;
proc sort data=event_sim nodup;
   by event_id;
run;

data junc_flag_sim;
   merge event2uniq (in=in1) event_sim (in=in2) event_gene_sim (in=in3);
   by event_id;
   if in2 then flag_junction_xscript_sim=1;
   else flag_junction_xscript_sim=0;
   if in2 then flag_junction_gene_sim=1;
   else flag_junction_gene_sim=0;
   if in1 then output;
   keep junction_id seq_name flag_junction_xscript_sim flag_junction_gene_sim;
run;

proc sort data=junc_flag_sim nodup;
  by seq_name junction_id;
run;

data cat_junc;
   set event.catalog_junctions_10000xs;
   keep sample_id event_id apn;
   rename event_id=seq_name;
run;

proc sort data=cat_junc;
   by seq_name sample_id;
proc transpose data=cat_junc out=cat_junc_sbys;
   by seq_name;
   id sample_id;
   var apn;
run;

proc sort data=cat_junc_sbys;
  by seq_name;
proc sort data=junc_flag_sim;
  by seq_name;
run;

data cat_junc_sbys2;
  merge junc_flag_sim (in=in1) cat_junc_sbys (in=in2);
  by seq_name;
  if in1 and in2;
  drop seq_name;
run;

/* Add in STAR junctions */

data star_junc;
  set event.star_junctions_10000xs;
  length junction_id $50.;
  if max_overhang lt 16 then delete;
  if strand=1 then strand_str="+"; else strand_str="-";
  donor_site=intron_start-1;
  acceptor_site=intron_stop;
  junction_id=catx(":",chr,donor_site,acceptor_site,strand_str);
  keep sample_id junction_id num_unique_mapped_reads;
  rename num_unique_mapped_reads=reads_mapped_star;
run;

proc sort data=star_junc;
   by junction_id sample_id;
proc transpose data=star_junc out=star_junc_sbys;
   by junction_id;
   id sample_id;
   var reads_mapped_star;
run;

data star_junc_sbys2;
  set star_junc_sbys;
  drop _NAME_;
  rename sample_01=star_sample_01 sample_02=star_sample_02 sample_03=star_sample_03
         sample_04=star_sample_04 sample_05=star_sample_05 sample_06=star_sample_06;
run;


proc sort data=star_junc_sbys2 nodup;
  by junction_id;
proc sort data=cat_junc_sbys2  nodup;
  by junction_id;
run;

data cat2star_junc_sbys;
  merge cat_junc_sbys2 (in=in1) star_junc_sbys2 (in=in2);
  by junction_id;
  if in1 and in2 then do;
      flag_in_catalog=1;
      flag_in_star=1;
      output; end;
  else if in1 then do;
      flag_in_catalog=1;
      flag_in_star=0;
      star_sample_01=0; star_sample_02=0; star_sample_03=0;
      star_sample_04=0; star_sample_05=0; star_sample_06=0;
      output; end;
  else do;
      flag_in_catalog=0;
      flag_in_star=1;
      sample_01=0; sample_02=0; sample_03=0;
      sample_04=0; sample_05=0; sample_06=0;
      flag_junction_xscript_sim=0;
      flag_junction_gene_sim=0;
      output; end;
run;


data junc_types;
   set event.event2star2pacbio_junc_table;
   keep junction_id junction_type;
run;

proc sort data=cat2star_junc_sbys;
   by junction_id;
proc sort data=junc_types;
   by junction_id;
run;

data cat2star_junc_sbys2;
  merge junc_types (in=in1) cat2star_junc_sbys (in=in2);
  by junction_id;
  if in2;
  if junction_type="" then junction_type="STAR only";
  drop _NAME_;
run;


/* Flag detected by sample, count, export */

* transpose so I can do this by sample;

proc sort data=cat2star_junc_sbys2;
   by junction_id junction_type;
proc transpose data=cat2star_junc_sbys2 out=cat2star_junc_ea;
   by junction_id junction_type;
   var sample_01 sample_02 sample_03 sample_04 sample_05 sample_06;
run;

proc transpose data=cat2star_junc_sbys2 out=cat2star_junc_star;
   by junction_id junction_type;
   var star_sample_01 star_sample_02 star_sample_03 star_sample_04 star_sample_05 star_sample_06;
run;

data cat2star_junc_ea2;
   set cat2star_junc_ea;
   rename COL1=event_apn _NAME_=sample_id;
run;

data cat2star_junc_star2;
   set cat2star_junc_star;
   length sample_id $9.;
   sample_id=compress(tranwrd(_NAME_,"star_",""));
   rename COL1=star_mapped_reads;
   drop _NAME_;
run;

data cat2star_junc_flags;
   set cat2star_junc_sbys2;
   keep junction_id flag_in_catalog flag_in_star flag_junction_xscript_sim flag_junction_gene_sim;
run;

data junc_annot junc_unannot junc_border junc_staronly;
   set cat2star_junc_sbys2;
   if junction_type="Annotated" then output junc_annot;
   else if junction_type="NIC/Unannotated" then output junc_unannot;
   else if junction_type="STAR only" then output junc_staronly;
   else if junction_type="border junction" then output junc_border;
   keep junction_id;
run;

proc sort data=cat2star_junc_ea2 nodup;
   by junction_id junction_type sample_id;
proc sort data=cat2star_junc_star2 nodup;
   by junction_id junction_type sample_id;
proc sort data=cat2star_junc_flags nodup;
   by junction_id;
proc sort data=junc_annot nodup;
   by junction_id;
proc sort data=junc_unannot nodup;
   by junction_id;
proc sort data=junc_border nodup;
   by junction_id;
proc sort data=junc_staronly nodup;
   by junction_id;
run;

proc means data=cat2star_junc_flags noprint;
   by junction_id;
   var flag_junction_xscript_sim flag_junction_gene_sim flag_in_catalog flag_in_star;
   output out=cat2star_junc_flags2(drop=_TYPE_ _FREQ_) max=;
run;


data junc_types;
  merge junc_annot (in=in1) junc_unannot (in=in2) junc_border (in=in3) junc_staronly (in=in4); 
   by junction_id;
    length junction_type $50.;
   if in1 then junction_type="Annotated";
   else if in2 then junction_type="NIC/Unannotated";
   else if in3 then junction_type="border junction";
   else if in4 then junction_type="STAR only";
run;

data cat2star_junc;
   merge cat2star_junc_ea2 (in=in1) cat2star_junc_star2 (in=in2);
   by junction_id junction_type sample_id;
   if in1 and in2;
run;

data cat2star_junc2;
   merge cat2star_junc (in=in1) cat2star_junc_flags2 (in=in2)  junc_types;
   by junction_id ;
   if in1 and in2 then output;
run;

/* Flag if detected and count junctions per sample */

data cat2star_junc3;
  set cat2star_junc2;
  if event_apn > 0 then flag_events_detected=1; else flag_events_Detected=0;
  if star_mapped_reads > 0 then flag_star_detected=1; else flag_star_detected=0;
run;

proc sort data=cat2star_junc3;
  by sample_id junction_type junction_id;
proc freq data=cat2star_junc3 noprint;
  by sample_id;
  tables junction_type*flag_in_catalog*flag_in_star*flag_events_detected*flag_star_detected*
         flag_junction_xscript_sim / out=junc_counts_by_sample;
    tables flag_junction_xscript_sim*flag_events_detected*flag_star_detected / out=junc_counts_by_sample_basic;
run;

proc freq data=cat2star_junc3 noprint;
  by sample_id;

run;

proc print data=junc_counts_by_sample;
  where sample_id="sample_01";
run;

proc print data=junc_counts_by_sample_basic;
  where sample_id="sample_01";
run;

/*
                                                           flag_
                          flag_in_ flag_in_ flag_events_   star_  flag_junction_
sample_id junction_type    catalog   star     detected   detected   xscript_sim  COUNT PERCENT

JUNCTIONS FROM SIMULATED TRANSCRIPTS
sample_01 Annotated           1        0          0          0           1          24  0.0392
sample_01 Annotated           1        0          1          0           1         425  0.6938
sample_01 Annotated           1        1          0          0           1          48  0.0784
sample_01 Annotated           1        1          0          1           1          20  0.0326
sample_01 Annotated           1        1          1          0           1         137  0.2236
sample_01 Annotated           1        1          1          1           1       57544 93.9371

48+20+24 from transcripts that aren't detected with events in this sample
24+48 of these are also missed by STAR in this sample

ALL OTHER
sample_01 Annotated           0        1          0          0           0         101  0.1649
sample_01 Annotated           0        1          0          1           0         155  0.2530
sample_01 Annotated           1        0          0          0           0          70  0.1143
sample_01 Annotated           1        0          1          0           0         837  1.3664
sample_01 Annotated           1        1          0          0           0          16  0.0261
sample_01 Annotated           1        1          0          1           0           6  0.0098
sample_01 Annotated           1        1          1          0           0          39  0.0637
sample_01 Annotated           1        1          1          1           0         130  0.2122
sample_01 NIC/Unannotated     0        1          0          0           0           7  0.0114
sample_01 NIC/Unannotated     0        1          0          1           0           4  0.0065
sample_01 NIC/Unannotated     1        0          0          0           0          50  0.0816
sample_01 NIC/Unannotated     1        0          1          0           0         358  0.5844
sample_01 NIC/Unannotated     1        1          1          0           0           2  0.0033
sample_01 STAR only           0        1          0          0           0         203  0.3314
sample_01 STAR only           0        1          0          1           0         133  0.2171
sample_01 border junction     1        0          0          0           0         139  0.2269
sample_01 border junction     1        0          1          0           0         810  1.3223

Simplified:

                                                         flag_
                     flag_junction_    flag_events_      star_
 Obs    sample_id      xscript_sim       detected      detected    COUNT    PERCENT

   1    sample_01           0                0             0         586     0.9566
   2    sample_01           0                0             1         298     0.4865
   3    sample_01           0                1             0        2046     3.3400
   4    sample_01           0                1             1         130     0.2122
   5    sample_01           1                0             0          72     0.1175
   6    sample_01           1                0             1          20     0.0326
   7    sample_01           1                1             0         562     0.9174
   8    sample_01           1                1             1       57544    93.9371

*/


proc export data=junc_counts_by_sample
     outfile="!MCLAB/event_analysis/analysis_output/catalog_vs_star_10000xscripts_junctions_by_sample.csv"
     dbms=csv replace;
run;



proc export data=junc_counts_by_sample_basic
     outfile="!MCLAB/event_analysis/analysis_output/catalog_vs_star_10000xscripts_junctions_by_sample_simple.csv"
     dbms=csv replace;
run;


