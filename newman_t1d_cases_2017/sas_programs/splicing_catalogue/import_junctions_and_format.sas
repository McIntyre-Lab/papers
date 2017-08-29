/***** IMPORT JUNCTIONS AND FORMAT *****/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';


/* Import junctions */

    data WORK.LOGICAL_JUNCTIONS_BED    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_junctions.bed'
delimiter='09'x MISSOVER DSD lrecl=32767 ;
       informat chr $4. ;
       informat totalstart best32. ;
       informat totalstop best32. ;
       informat event_id $79. ;
       informat score best32. ;
       informat strand $1. ;
       informat totalstart2 best32. ;
       informat totalstop2 best32. ;
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
       format totalstart2 best12. ;
       format totalstop2 best12. ;
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
               totalstart2
               totalstop2
               color $
               blocks
               block_sizes $
               block_starts $
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;


/* Format junctions: need donor and acceptor exons, and get all start and stop positions */

data junctions_formatted;
   set logical_junctions_bed;
   length donor_exon $39.;
   length acceptor_exon $39.;
   length event_type $16.;
   donor_exon=scan(event_id,1,'|');
   acceptor_exon=scan(event_id,2,'|');
   donor_size=scan(block_sizes,1,',')+0;
   acceptor_size=scan(block_sizes,2,',')+0;
   donor_start=totalstart;
   donor_stop=totalstart+donor_size;126
   acceptor_start=totalstop-acceptor_size;
   acceptor_stop=totalstop;
   event_type='exon_junction';
   drop totalstart totalstop totalstart2 totalstop2 score color blocks block_sizes block_starts;
run;


/* Save as permenant */

data splice.logical_junctions;
   set junctions_formatted;
run;

/* remove unwanted datasets */
proc datasets noprint;
   delete junctions_formatted logical_junctions_bed;
run;
quit;

