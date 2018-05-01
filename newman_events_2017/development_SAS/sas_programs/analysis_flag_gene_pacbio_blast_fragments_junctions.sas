ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Process PacBio BLAST results:
- Drop hits where the hit length is <90% of the event's length
- Flag if 90%, 95%, 99%, 100%
- Flag if event is going to at least the correct gene
- Flag if event is going to at least one of its originating transcripts
- Flag if event if only going to another gene and not its own

Remove features not expressed
Remove redundant features
Remove features from genes with multigene components
Remove hits from genes not expressed
Check feature is (at least) going to the correct gene, and flag if not
*/

data blast_w_feat_length;
   set event.pacbio_blast_hits_to_parse;
run;

data drop_short_events;
  set blast_w_feat_length;
  if length < 0.9 * feature_length then delete;
  if length >= 0.9 * feature_length then flag_length_90perc=1; else flag_length_90perc=0;
  if length >= 0.95 * feature_length then flag_length_95perc=1; else flag_length_95perc=0;
  if length >= 0.99 * feature_length then flag_length_99perc=1; else flag_length_99perc=0;
run;

proc freq data=drop_short_events;
  tables flag_length_90perc flag_length_95perc flag_length_99perc;
run;

/*

       flag_length_                             Cumulative    Cumulative
             90perc    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  1      140303      100.00        140303       100.00


       flag_length_                             Cumulative    Cumulative
             95perc    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0         675        0.48           675         0.48
                  1      139628       99.52        140303       100.00


       flag_length_                             Cumulative    Cumulative
             99perc    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        2299        1.64          2299         1.64
                  1      138004       98.36        140303       100.00


I can probably drop any that has a length < 99% of feature length, as I'm not losing too much.
Leave flags for now and decide later

*/

/* Merge in RefSeq ID for PacBio transcript */

data pb2xs;
   set event.pacbio2refseq_id_nomulti;
   keep pacbio_id transcript_id;
run;

proc sort data=pb2xs nodup;
   by pacbio_id transcript_id;
proc sort data=drop_short_events;
   by pacbio_id;
run;

data blast_w_pb2xs;
  merge pb2xs (in=in1) drop_short_events (in=in2);
  by pacbio_id;
  if in2;
run;

/* Flag if event gene is the same as PacBio gene */

data flag_gene;
  set blast_w_pb2xs;
  if count(refseq_gene_id,strip(gene_id)) > 0 then flag_same_gene=1;
  else flag_same_gene=0;
run;

proc freq data=flag_gene;
  tables flag_same_gene;
run;


/*
                                            Cumulative    Cumulative
 flag_same_gene    Frequency     Percent     Frequency      Percent
 -------------------------------------------------------------------
              0        1172        0.84          1172         0.84
              1      139131       99.16        140303       100.00


Try at the event-level. I want to know if the event is going to the correct gene at all
*/

proc sort data=flag_Gene;
   by feature_id;
proc means data=flag_gene noprint;
   by feature_id;
   var flag_same_gene; 
   output out=sum_flag_gene max=;
run;

proc freq data=sum_flag_gene;
  table flag_same_gene;
run;

/*
                                           Cumulative    Cumulative
flag_same_gene    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0         121        0.18           121         0.18
             1       65659       99.82         65780       100.00

Okay, there are some hits that are not going to the correct gene. I am going to put these into their own dataset for looking at later

*/

data same_gene diff_gene;
   set flag_gene;
   if flag_same_gene=1 then output same_gene;
   else output diff_gene;
run;

/* one last check: how many are hits to known transcripts? */

data fusxs_check;
   set diff_gene;
   if pacbio_status = "Known";
run; *643 of 1172 hits;


/* Make permenant */
data event.pacbio_hit_wrong_gene;
  set diff_gene;
run;

data event.pacbio_hit_correct_gene;
  set same_gene;
run;


