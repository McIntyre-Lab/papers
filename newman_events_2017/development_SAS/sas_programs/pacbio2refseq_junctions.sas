/* Import Pacbio transcript junctions and figure out which are detected, corresponding to the catalog

I have pulled the donor and acceptor coordinates per transcript, so I need to match this to feature1_stop and feature2_start
*/

libname event '!MCLAB/event_analysis/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';

data events;
   set evspl.splicing_events_annot_refseq;
   if num_transcripts > 0 then flag_junction_annotated=1;
   keep event_id gene_id chr strand transcript_id
        feature1_start feature1_stop feature2_start feature2_stop
        flag_junction_annotated flag_intron_retention flag_exonskip
        flag_alt_donor flag_alt_acceptor;
run;

data on_events;
  set refseq.rfsq_flag_splicing_on;
  keep event_id flag_splicing_on;
run;

data counts;
   set refseq.rfsq_counts_by_splicing;
run;

proc sort data=counts;
  by event_id;
proc means data=counts noprint;
  by event_id;
  var apn;
  output out=mean_counts mean= / autoname;
run;

*import pacbio transcript junctions;

proc import datafile="!MCLAB/event_analysis/analysis_output/mm10_conesa_pacbio_junctions.csv" out=pb_junctions
     dbms=csv replace; guessingrows=145172;
run;

*format junctions;

data pb_junctions2;
   set pb_junctions;
   rename chrom=pb_event_id start=coord stop=pb_transcript_id strand=pb_gene_id;
run;

*Cat event Ids;
data junc2id;
  set pb_junctions2;
  keep coord pb_event_id;
run;

proc sort data=junc2id nodup;
  by coord;
proc freq data=junc2id noprint;
  tables coord / out=coord_count;
proc sort data=coord_count;
  by descending count;
proc print data=coord_count (obs = 1);
run; *26 transcripts per coord;

data cat_id; 
  array xs[26] $ 31.;
  retain xs1-xs26;
  set junc2id;
  by coord;
  if first.coord then do;
     call missing(of xs1-xs26);
     records = 0;
  end;
  records + 1;
  xs[records]=pb_event_id;
  if last.coord then output;
run;

data cat_id2;
  set cat_id;
  length cat_event_id $ 850.;
         cat_event_id= catx("|", OF xs1-xs26);
  keep coord cat_event_id records;
  rename records= num_events;
  rename cat_event_id=pb_event_id;
  run;


* Cat transcripts;
data junc2xs;
  set pb_junctions2;
  keep coord pb_transcript_id;
run;

proc sort data=junc2xs nodup;
  by coord;
proc freq data=junc2xs noprint;
  tables coord / out=coord_count;
proc sort data=coord_count;
  by descending count;
proc print data=coord_count (obs = 1);
run; *33 transcripts per coord;

data cat_xs; 
  array xs[33] $ 15.;
  retain xs1-xs33;
  set junc2xs;
  by coord;
  if first.coord then do;
     call missing(of xs1-xs33);
     records = 0;
  end;
  records + 1;
  xs[records]=pb_transcript_id;
  if last.coord then output;
run;

data cat_xs2;
  set cat_xs;
  length cat_xs_id $ 550.;
         cat_xs_id= catx("|", OF xs1-xs33);
  keep coord cat_xs_id records;
  rename records= num_transcripts;
  rename cat_xs_id=pb_transcript_id;
  run;


* Cat genes;
data junc2gene;
  set pb_junctions2;
  keep coord pb_gene_id;
run;

proc sort data=junc2gene nodup;
  by coord;
proc freq data=junc2gene noprint;
  tables coord / out=coord_count;
proc sort data=coord_count;
  by descending count;
proc print data=coord_count (obs = 1);
run; *6 gene per coord;

data cat_gene; 
  array xs[6] $ 15.;
  retain xs1-xs6;
  set junc2gene;
  by coord;
  if first.coord then do;
     call missing(of xs1-xs6);
     records = 0;
  end;
  records + 1;
  xs[records]=pb_gene_id;
  if last.coord then output;
run;

data cat_gene2;
  set cat_gene;
  length cat_gene_id $ 90.;
         cat_gene_id= catx("|", OF xs1-xs6);
  keep coord cat_gene_id records;
  rename records= num_genes;
  rename cat_gene_id=pb_gene_id;
  run;

proc sort data=cat_id2;
  by coord;
proc sort data=cat_xs2;
  by coord;
proc sort data=cat_gene2;
  by coord;
run;

data pb_junctions3;
  merge cat_id2 cat_xs2 cat_gene2;
  by coord;
run;


data pb_junctions4;
  length chr $5.;
  format feature1_stop best32. ;
  format feature2_start best32. ;
  length strand $1.;
  set pb_junctions3;
  chr=scan(coord,1,":");
  feature1_stop=scan(coord,2,":");
  feature2_start=scan(coord,3,":") - 1;
  strand=scan(coord,4,":");
run;

proc sort data=pb_junctions4;
  by chr feature1_stop feature2_start strand;
proc sort data=events;
  by chr feature1_stop feature2_start strand;
run;

data pb2refseq noref;
   merge events (in=in1) pb_junctions4 (in=in2);
   by chr feature1_stop feature2_start strand;
   if in1 then output pb2refseq;
   else if in2 then output noref;
run;

/* 3691 PB junctions have no matching RefSeq junction */

* How many PB junctions with a Refseq catalog match are
  annotated junctions
  unannotated junctions
  IR events
  exonskipping
  alt donor
  alt acceptor
  alt donor and acceptor ;

data pb2refseq_only;
  set pb2refseq;
  if pb_event_id ne '';
run; *69045 catalog events have a PB match;

proc freq data=pb2refseq_only;
  tables flag_junction_annotated*flag_intron_retention;
run;

*68393 annotated junctions;
*652 unannotated junctions;
*0 IR events, expect this;

proc freq data=pb2refseq_only;
  where flag_intron_retention=0;
  tables flag_exonskip flag_alt_donor*flag_alt_acceptor;
run;

*5995 exonskipping;
*1512 alt donor;
*1614 alt acceptor;
*514 alt donor+acceptor;

/* Merge in on flags and mean counts. Let's set APN thresholds of 5, 10 and 20, and count the number of PB junctions */

proc sort data=pb2refseq;
  by event_id;
proc sort data=on_events;
  by event_id;
proc sort data=mean_counts;
  by event_id;
run;

data pb2refseq_w_flags;
   merge pb2refseq (in=in1) on_events (in=in2) mean_counts (in=in3);
   by event_id;
   if apn_mean < 5 then flag_apn_ge5=0; else flag_apn_ge5=1;
   if apn_mean < 10 then flag_apn_ge10=0; else flag_apn_ge10=1;
   if apn_mean < 20 then flag_apn_ge20=0; else flag_apn_ge20=1;
   if pb_event_id='' then flag_pb_junction=0; else flag_pb_junction=1;
run;

/* Check: how many PB transcripts are detected at APN>0, 5, 10, 20? */

proc freq data=pb2refseq_w_flags;
   tables flag_splicing_on*flag_pb_junction
          flag_apn_ge5*flag_pb_junction
          flag_apn_ge10*flag_pb_junction
          flag_apn_ge20*flag_pb_junction;
run;

/* 64247 of 69045 catalog events with a PB junction are detected at APN>0 
 44709 of 69045 catalog events with a PB junction are detected at APN>5
 34245 of 69045 catalog events with a PB junction are detected at APN>10 
 23604 of 69045 catalog events with a PB junction are detected at APN>20
*/

/* Of the 44709 events with PB junctions detected at APN>0/5/10/20, how many are not annotated? (will be a max of 652)*/

proc freq data=pb2refseq_w_flags;
   where flag_pb_junction=1;
   tables flag_splicing_on*flag_junction_annotated
          flag_apn_ge5*flag_junction_annotated
          flag_apn_ge10*flag_junction_annotated
          flag_apn_ge20*flag_junction_annotated;
run;

/* At APN>0, there are 327 of 652 unannotated junctions with a PB counterpart
   At APN>5, there are 48 of 652 unannotated junctions with a PB counterpart
   At APN>10, there are 20 of 652 unannotated junctions with a PB counterpart
   At APN>20, there are 11 of 652 unannotated junctions with a PB counterpart

   Are transcripts with unannotated events lowly expressed? */
