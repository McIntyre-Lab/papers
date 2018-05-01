/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the 59-gene simulation, I want to count the number of times
   we see a junction in the catalog, in STAR, the union, and the intersect.
   Do by sample, and then for those junctions seen in all samples */

data annot;
  set evspl.splicing_events_annot_refseq;
  length junction_id $50.;
  where num_transcripts > 0 or flag_junction_annotated=1;
  junction_id=catx(":",chr,feature1_stop,feature2_start,strand);
  keep junction_id transcript_id;
run;

data annot2;
  set annot;
  length transcript_id2 $20.;
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "" );
    transcript_id2=scan(transcript_id,i,"|");
    output;
    end;
  drop transcript_id;
  rename transcript_id2=transcript_id;
run;

/* Get list of transcripts simulated */

data xs2gene;
  set event.feature2xs2gene;
  keep transcript_id gene_id;
run;

data genelist_1;
  format gene_id $11.;
  input gene_id $;
  datalines;
  101214
  106581
  108011
  114565
  114863
  11864
  12419
  12830
  142980
  16150
  16468
  170787
  17350
  20509
  217378
  22029
  22276
  225608
  22654
  228960
  231386
  232337
  234366
  234725
  269120
  27416
  28240
  56812
  58175
  66433
  66855
  66940
  67669
  69930
  70380
  70382
  70796
  71893
  72205
  72344
  72999
  73750
  74080
  74140
  74190
  76577
  77613
  77733
  94281
  ;
run;

data genelist_2;
  set event.genes_w_nic_junction_10genes;
run;

proc sort data=genelist_1;
  by gene_id;
proc sort data=genelist_2;
  by gene_id;
run;

data genelist_for_sim;
  merge genelist_1 (in=in1) genelist_2 (in=in2);
  by gene_id;
  if in2 then flag_has_nic=1; else flag_has_nic=0;
run;

proc sort data=genelist_for_sim;
  by gene_id;
proc sort data=xs2gene nodup;
  by gene_id;
run;

data xslist_sim;
  merge genelist_for_sim (in=in1) xs2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

proc sort data=annot2 nodup;
   by transcript_id  junction_id;
proc sort data=xslist_sim;
   by transcript_id;
run;

data annot_sim;
  merge annot2 (in=in1) xslist_sim (in=in2);
  by transcript_id;
  if in1 and in2;
  keep junction_id;
run;
proc sort data=annot_sim nodup;
  by junction_id;
run;
*1075 junctions from these 467 transcripts/59 genes;

/* How many of these junctions are detected with STAR and via the catalog? */

data cat_junc;
   set event.catalog_junctions_10genes;
   if apn>0 then flag_detected=1; else flag_detected=0;
   keep sample_id event_id flag_detected;
   rename event_id=seq_name;
run;

proc sort data=cat_junc;
   by seq_name sample_id;
proc transpose data=cat_junc out=cat_junc_sbys;
   by seq_name;
   id sample_id;
   var flag_detected;
run;

data cat_seq2id;
  set eventloc.unique_junction2event_mm10;
  keep seq_name junction_id;
run;

proc sort data=cat_junc_sbys;
  by seq_name;
proc sort data=cat_seq2id nodup;
  by seq_name;
run;

data cat_junc_id;
  merge cat_junc_sbys (in=in1) cat_seq2id (in=in2);
  by seq_name;
  if in1 and in2;
  drop seq_name;
run;

proc sort data=cat_junc_id;
  by junction_id;
proc sort data=annot_sim;
  by junction_id;
run;

data cat_junc_id_sim;
  merge annot_sim (in=in1) cat_junc_id (in=in2);
  by junction_id;
  if in1 and in2 then output;
  else if in1 then do;
    sample_01=0; sample_02=0; sample_03=0;
    sample_04=0; sample_05=0; sample_06=0;
    output;
    end;
  drop _NAME_;
  rename sample_01=flag_sample1_events_dtct
         sample_02=flag_sample2_events_dtct
         sample_03=flag_sample3_events_dtct
         sample_04=flag_sample4_events_dtct
         sample_05=flag_sample5_events_dtct
         sample_06=flag_sample6_events_dtct;
run;

/* STAR junctions */


data star_junc;
  set event.star_junctions_10genes;
  length junction_id $50.;
  if max_overhang lt 16 then delete;
  if strand=1 then strand_str="+"; else strand_str="-";
  donor_site=intron_start-1;
  acceptor_site=intron_stop;
  junction_id=catx(":",chr,donor_site,acceptor_site,strand_str);
  if num_unique_mapped_reads > 0 then flag_detected=1;
  else flag_detected=0;
  keep sample_id junction_id flag_detected;
run;

proc sort data=star_junc;
   by junction_id sample_id;
proc transpose data=star_junc out=star_junc_sbys;
   by junction_id;
   id sample_id;
   var flag_detected;
run;

data star_junc_sbys2;
  set star_junc_sbys;
  drop _NAME_;
  rename sample_01=flag_sample1_star_dtct
         sample_02=flag_sample2_star_dtct
         sample_03=flag_sample3_star_dtct
         sample_04=flag_sample4_star_dtct
         sample_05=flag_sample5_star_dtct
         sample_06=flag_sample6_star_dtct;
run;

data star_junc_sbys3;
  set star_junc_sbys2;
  array change _numeric_;
  do over change;
   if change=. then change=0;
  end;
run;



proc sort data=star_junc_sbys3 nodup;
  by junction_id;
proc sort data=cat_junc_id_sim nodup;
  by junction_id;
run;

data cat2star_sim;
  merge cat_junc_id_sim (in=in1) star_junc_sbys3 (in=in2);
  by junction_id;
  if in1 and in2 then output;
  else if in1 then do;
       flag_sample1_star_dtct=0;
       flag_sample2_star_dtct=0;
       flag_sample3_star_dtct=0;
       flag_sample4_star_dtct=0;
       flag_sample5_star_dtct=0;
       flag_sample6_star_dtct=0;
       output; end;
run;


/* Flag if detected in all samples */

data cat2star_sim_flag;
  set cat2star_sim;
  if sum(flag_sample1_events_dtct,flag_sample2_events_dtct,flag_sample3_events_dtct,
         flag_sample4_events_dtct,flag_sample5_events_dtct,flag_sample6_events_dtct) = 6
     then flag_all_events_dtct=1; else flag_all_events_dtct=0;

  if sum(flag_sample1_star_dtct,flag_sample2_star_dtct,flag_sample3_star_dtct,
         flag_sample4_star_dtct,flag_sample5_star_dtct,flag_sample6_star_dtct) = 6
     then flag_all_star_dtct=1; else flag_all_star_dtct=0;

  if sum(flag_sample1_events_dtct,flag_sample2_events_dtct,flag_sample3_events_dtct,
         flag_sample4_events_dtct,flag_sample5_events_dtct,flag_sample6_events_dtct) > 0
     then flag_any_events_dtct=1; else flag_any_events_dtct=0;

  if sum(flag_sample1_star_dtct,flag_sample2_star_dtct,flag_sample3_star_dtct,
         flag_sample4_star_dtct,flag_sample5_star_dtct,flag_sample6_star_dtct) > 0
     then flag_any_star_dtct=1; else flag_any_star_dtct=0;
run;

/* Count */

proc freq data=cat2star_sim_flag noprint;
  tables flag_sample1_events_dtct*flag_sample1_star_dtct / out=sample1_events_v_star;
  tables flag_sample2_events_dtct*flag_sample2_star_dtct / out=sample2_events_v_star;
  tables flag_sample3_events_dtct*flag_sample3_star_dtct / out=sample3_events_v_star;
  tables flag_sample4_events_dtct*flag_sample4_star_dtct / out=sample4_events_v_star;
  tables flag_sample5_events_dtct*flag_sample5_star_dtct / out=sample5_events_v_star;
  tables flag_sample6_events_dtct*flag_sample6_star_dtct / out=sample6_events_v_star;
  tables flag_all_events_dtct*flag_all_star_dtct / out=all_events_v_star;
  tables flag_any_events_dtct*flag_any_star_dtct / out=any_events_v_star;
run;

proc print data=sample1_events_v_star;
run;
proc print data=sample2_events_v_star;
run;
proc print data=sample3_events_v_star;
run;
proc print data=sample4_events_v_star;
run;
proc print data=sample5_events_v_star;
run;
proc print data=sample6_events_v_star;
run;
proc print data=all_events_v_star;
run;
proc print data=any_events_v_star;
run;

/* Tables 

SAMPLE 1:
 flag_sample1_    flag_sample1_
  events_dtct       star_dtct      COUNT

       0                0             12
       0                1              1
       1                0             18
       1                1           1044


SAMPLE 2:
flag_sample2_    flag_sample2_
 events_dtct       star_dtct      COUNT

      0                0             12
      0                1              1
      1                0             18
      1                1           1044



SAMPLE 3:
 flag_sample3_    flag_sample3_
  events_dtct       star_dtct      COUNT

       0                0             12
       1                0             18
       1                1           1045

SAMPLE 4:
flag_sample4_    flag_sample4_
 events_dtct       star_dtct      COUNT

      0                0             12
      1                0             18
      1                1           1045


SAMPLE 5:
flag_sample5_    flag_sample5_
 events_dtct       star_dtct      COUNT

      0                0             13
      1                0             17
      1                1           1045


SAMPLE 6:
flag_sample6_    flag_sample6_
 events_dtct       star_dtct      COUNT

      0                0             12
      1                0             18
      1                1           1045

ALL SAMPLES:
flag_all_
 events_     flag_all_
   dtct      star_dtct    COUNT

    0            0           13
    0            1            1
    1            0           19
    1            1         1042


ANY SAMPLE:
flag_any_
 events_     flag_any_
   dtct      star_dtct    COUNT

    0            0           12
    1            0           17
    1            1         1046

*/


