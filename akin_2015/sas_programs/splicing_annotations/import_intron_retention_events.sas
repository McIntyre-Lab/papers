libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

/* Import IR events */


     data WORK.INTRON_RETENTION_EVENTS    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_intron_retention.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat gene_id $36. ;
        informat event_id $79. ;
        informat chr $4. ;
        informat strand $1. ;
        informat intron_position best32. ;
        informat exon_id $39. ;
        informat exon_cat $1283. ;
        informat flag_lastexon best32. ;
        format gene_id $36. ;
        format event_id $79. ;
        format chr $4. ;
        format strand $1. ;
        format intron_position best12. ;
        format exon_id $39. ;
        format exon_cat $1283. ;
        format flag_lastexon best12. ;
     input
                 gene_id $
                 event_id $
                 chr $
                 strand $
                 intron_position
                 exon_id $
                 exon_cat $
                 flag_lastexon
     ;
     if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
     run;



     data WORK.INTRON_RETENTION_EVENTS_BED    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_intron_retention.bed' delimiter='09'x MISSOVER DSD lrecl=32767 ;
        informat chr $4. ;
       informat totalstart best32. ;
       informat totalstop best32. ;
       informat event_id $79. ;
       informat score best32. ;
       informat strand $1. ;
       informat color $8. ;
       informat blocks best32. ;
       informat block_sizes $6. ;
       informat block_starts $9. ;
       format chr $4. ;
       format totalstart best12. ;
       format totalstop best12. ;
       format event_id $79. ;
       format score best12. ;
       format strand $1. ;
       format color $8. ;
       format blocks best12. ;
       format block_sizes $6. ;
       format block_starts $9. ;
        input
               chr $
               totalstart
               totalstop
               event_id $
               score
               strand $
               color $
               blocks
               block_sizes $
               block_starts $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;



/* Format datasets */


data ir_events_formatted;
   set INTRON_RETENTION_EVENTS_BED;
   length donor_exon $39.;
   length acceptor_exon $39.;
   length event_type $16.;
   if strand='+' then do;
       donor_exon=scan(event_id,1,'|');
       acceptor_exon=scan(event_id,2,'|');
       event_size=block_sizes+0;
       donor_size=event_size-37;
       acceptor_size=36;
       donor_start=totalStart;
       donor_stop=totalStart+donor_size;
       acceptor_start=totalStop-acceptor_size;
       acceptor_stop=totalStop;
       end;
   else if strand='-' then do;
       donor_exon=scan(event_id,2,'|');
       acceptor_exon=scan(event_id,1,'|');
       event_size=block_sizes+0;
       donor_size=36;
       acceptor_size=event_size-37;
       donor_start=totalStart;
       donor_stop=totalStart+donor_size;
       acceptor_start=totalStop-acceptor_size;
       acceptor_stop=totalStop;
       end;
    event_type='intron_retention';
    flag_intron_retention=1;
   drop totalstart totalstop score color blocks block_sizes block_starts;
run;

/* Make permenant */

data splice.intron_retention_info;
   set intron_retention_events;
run;

data splice.intron_retention_events;
   set ir_events_formatted;
run;


/* remove temp datasets */

proc datasets noprint;
  delete INTRON_RETENTION_EVENTS_BED intron_retention_events ir_events_formatted;
run;
quit;

