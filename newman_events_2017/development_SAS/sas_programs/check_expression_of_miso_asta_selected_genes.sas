/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Check: the genes selected for testing MISO: are they expresssed? Are the events expressed? */

data genes2test;
  input gene_id $11.;
  datalines;
  21419
  217378
  234725
  24058
  435802
  666528
  70375
  71893
  74140
  74190
;
run;

/* Get junctions */

data events2check;
  set evspl.splicing_events_annot_refseq;
  where gene_id in ('21419','217378','234725','24058','435802',
                    '666528','70375','71893','74140','74190');
  if num_transcripts > 0 then flag_junction_annotated=1;
  keep event_id gene_id transcript_id num_transcripts flag_junction_annotated flag_intron_retention
       flag_exonskip flag_alt_donor flag_alt_acceptor;
run;

data events_on;
  set eventloc.unique_junction2event_mm10;
  drop seq_name;
run;

proc sort data=events2check;
   by event_id;
proc sort data=events_on;
   by event_id;
run;

data events2check_2;
  merge events2check (in=in1) events_on (in=in2);
  by event_id;
  if in1 and in2;
run;

proc sort data=events2check_2;
  by gene_id;
proc freq data=events2check_2 noprint;
  by gene_id;
  tables flag_junction_annotated*flag_intron_retention*flag_exonskip*
         flag_alt_donor*flag_alt_Acceptor*flag_event_nsc_on / out=junc_by_gene_count;
run;

/* looks like the genes I picked aren't really expressed. I am going to re-run the gene selection code,
  limiting to only those genes with at least 5 annotated junctions detected */
