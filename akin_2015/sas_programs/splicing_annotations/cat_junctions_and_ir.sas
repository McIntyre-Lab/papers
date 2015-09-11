/***** Concatenate intron retention events *******/

libname splice '!MCLAB/junction_annotations/sas_data';


/* Format intron retention events so they can be appended to junctions */

/* Junctions columns: chr, event_id, strand, donor_exon, acceptor_exon, event_type, donor_size, acceptor_size, donor_start
   donor_stop, acceptor_start, acceptor_stop, num_transcripts, transcript_id, flag_junctions_annotated flag_exon_skip, num_skipped_exons */

/* IR columns: chr, event_id, strand, donor_exon, acceptor_exon, event_type, event_size, donor_size, acceptor_size, donor_start, donor_stop
   acceptor_start, acceptor_stop flag_intron_retention */

data intron_retention_events;
    set splice.intron_retention_events;
    num_transcripts=0;
    length transcript_id $2970.;
    transcript_id=' ';
    flag_junction_annotated=0;
    flag_exonskip=0;
    num_skipped_exons=0;
    drop event_size;
run;

data junctions_and_ir_events;
   set splice.junctions_w_exonskip (in=in1) intron_retention_events;
   if in1 then flag_intron_retention=0;
run;

/* Make permenant */

data splice.junction_and_ir_events;
    set junctions_and_ir_events;
run;

/* Delete temp datasets */

proc datasets noprint;
  delete intron_retention_events junctions_and_ir_events;
run;
quit;

