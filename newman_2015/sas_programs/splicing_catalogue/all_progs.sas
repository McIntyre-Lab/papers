/***** IMPORT JUNCTIONS AND FORMAT *****/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';


/* Import junctions */

    data WORK.LOGICAL_JUNCTIONS_BED    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_junctions_176bp.bed'
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
   donor_stop=totalstart+donor_size;
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




/***** IMPORT EXONS AND FORMAT *****/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

/* Import exons */

    data WORK.EXONS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_exons.csv'
delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat chrom $4. ;
       informat start best32. ;
       informat stop best32. ;
       informat strand $1. ;
       informat exon_id $39. ;
       informat transcript_id $1522. ;
       informat gene_id $36. ;
       format chrom $4. ;
       format start best12. ;
       format stop best12. ;
       format strand $1. ;
       format exon_id $39. ;
       format transcript_id $1522. ;
       format gene_id $36. ;
    input
                chrom $
                start
                stop
                strand $
                exon_id $
                transcript_id $
                gene_id $
    ;
  if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
  run;
   
/* Make permenant set - can use this later for things */


data av.hg19_aceview_exons;
   set exons;
   exon_length=stop-start;
run;

/* format exons and calculate length */

data exons_formatted;
   set exons;
   exon_length=stop-start;
   if exon_length lt 37 then flag_short_exon=1;
   else flag_short_exon=0;
   drop transcript_id;
run;


/* Group exons */
/* Sort exons by start and end
   Calculate length
   Group exons (see IR script for logic)
   Flag longest
   If longest, reference donor/acceptor
   If not longest, alternative donor/acceptor */


proc sort data=exons_formatted;
   by gene_id chrom start stop;
run;

/* test this with PINK1 and some other genes */


data pink1_exons;
  set exons_formatted;
  if gene_id='PINK1' or gene_id='PTPN22' or gene_id='GCH1' or gene_id='UBASH3A' or gene_id='tutee';
run;


data exon_grouping;
   retain exon_group region_start region_stop;
   set pink1_exons;
   by gene_id;
   if first.gene_id then do;
       exon_group=1;
       region_start=start;
       region_stop=stop;
       end;
   else do;
       if start le region_stop then do; *check to see if exon overlaps with previous exon;
           * if exon overlaps previous see which has the bigger coordinates;
           if stop ge region_stop then do; *expand exon region;
                region_start=region_start;
                region_stop=stop;
                end;
           else do; * exon is inside than current region's coordinates;
                region_start=region_start;
                region_stop=region_stop;
                end;
           exon_group=exon_group;
           end;
       else do; * exon does not overlap;
           exon_group=exon_group+1;
           region_start=start;
           region_stop=stop;
           end;
       end;
run;

/* test works, run on full exon set */

data exon_grouping;
   retain exon_group region_start region_stop;
   set exons_formatted;
   by gene_id;
   if first.gene_id then do;
       exon_group=1;
       region_start=start;
       region_stop=stop;
       end;
   else do;
       if start le region_stop then do; *check to see if exon overlaps with previous exon;
           * if exon overlaps previous see which has the bigger coordinates;
           if stop ge region_stop then do; *expand exon region;
                region_start=region_start;
                region_stop=stop;
                end;
           else do; * exon is inside the current region's coordinates;
                region_start=region_start;
                region_stop=region_stop;
                end;
           exon_group=exon_group;
           end;
       else do; * exon does not overlap with exon region in memory, then begin a new exon region;
           exon_group=exon_group+1;
           region_start=start;
           region_stop=stop;
           end;
       end;
run;


/* Seems to work! */

/* Flag longest exon in group as reference exon */
/* If two are the longest, then make the first of the two the reference exon */

proc sort data=exon_grouping;
   by gene_id exon_group descending exon_length start;
run;

data set_exon_reference;
   set exon_grouping;
   by gene_id exon_group;
   if first.exon_group then flag_reference_exon=1;
   else flag_reference_exon=0;
run;

/* Need to get references for each group */

proc sort data=set_exon_reference;
   by gene_id exon_group region_start region_stop;
run;

data exon_references;
   set set_exon_reference;
   by gene_id exon_group;
   if flag_reference_exon=1;
   keep gene_id exon_group start stop;
   rename start=ref_start stop=ref_stop;
run;

proc sort data=exon_references;
   by gene_id exon_group;
run;

proc sort data=set_exon_reference;
  by gene_id exon_group;
run;

data exon_groups_w_refs oops1 oops2;
   merge exon_references (in=in1) set_exon_reference (in=in2);
   by gene_id exon_group;
   if in1 and in2 then output exon_groups_w_refs;
   else if in1 then output oops1;
   else output oops2;
run;

/* Add in alternative donor/acceptor flags */

data exon_flag_alt_donor_acceptor;
   set exon_groups_w_refs;
   if start=ref_start then flag_alt_acceptor=0;
   else flag_alt_acceptor=1;
   if stop=ref_stop then flag_alt_donor=0;
   else flag_alt_donor=1;
run;

/* Checks to make sure that no reference exons have an alt donor or acceptor */
proc freq data=exon_flag_alt_donor_acceptor;
   tables flag_reference_exon*flag_alt_donor;
run;


/*

  flag_reference_exon
            flag_alt_donor

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 100277 | 228823 | 329100
           |  14.78 |  33.72 |  48.49
           |  30.47 |  69.53 |
           |  22.29 | 100.00 |
  ---------+--------+--------+
         1 | 349564 |      0 | 349564
           |  51.51 |   0.00 |  51.51
           | 100.00 |   0.00 |
           |  77.71 |   0.00 |
  ---------+--------+--------+
  Total      449841   228823   678664
              66.28    33.72   100.00

*/
proc freq data=exon_flag_alt_donor_acceptor;
   tables flag_reference_exon*flag_alt_acceptor;
run;

/*
 flag_reference_exon
           flag_alt_acceptor

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  99105 | 229995 | 329100
          |  14.60 |  33.89 |  48.49
          |  30.11 |  69.89 |
          |  22.09 | 100.00 |
 ---------+--------+--------+
        1 | 349564 |      0 | 349564
          |  51.51 |   0.00 |  51.51
          | 100.00 |   0.00 |
          |  77.91 |   0.00 |
 ---------+--------+--------+
 Total      448669   229995   678664
             66.11    33.89   100.00

*/

/* looks good! */

/* Make donor exon set */

data donor_exons;
   set exon_flag_alt_donor_acceptor;
   keep exon_group start stop exon_id gene_id flag_short_exon flag_alt_donor;
   rename exon_group=donor_group
          start=donor_exon_start
          stop=donor_exon_stop
          exon_id=donor_exon
          gene_id=donor_gene
          flag_short_exon=flag_short_donor;
run;

/* Make acceptor exon set */

data acceptor_exons;
   set exon_flag_alt_donor_acceptor;
   keep exon_group start stop exon_id gene_id flag_short_exon flag_alt_acceptor;
   rename exon_group=acceptor_group
          start=acceptor_exon_start
          stop=acceptor_exon_stop
          exon_id=acceptor_exon
          gene_id=acceptor_gene
          flag_short_exon=flag_short_acceptor;
run;

/* Make permenant */

data splice.donor_exons;
   set donor_exons;
run;

data splice.acceptor_exons;
  set acceptor_exons;
run;

data splice.exon_list;
    set exon_flag_alt_donor_acceptor;
run;

/* remove temp datasets */

proc datasets noprint;
  delete EXONS exons_formatted pink1_exons exon_grouping
  set_exon_reference  exon_flag_alt_donor_acceptor
  donor_exons acceptor_exons;
run;
quit;

/******** IMPORT TRANSCRIPT-ANNOTATED JUNCTIONS  **********/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

/* Import transcript-annotated junctions */

    data WORK.TRANSCRIPT_JUNCTIONS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_transcript_junctions.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat junction_id $79. ;
       informat junc_coords $24. ;
       informat transcript_id $53. ;
       informat gene_id $36. ;
       format junction_id $79. ;
       format junc_coords $24. ;
       format transcript_id $53. ;
       format gene_id $36. ;
    input
                junction_id $
                junc_coords $
                transcript_id $
                gene_id $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* Format transcript junctions */

data xscript_junc_formatted;
   set transcript_junctions;
   drop junc_coords; * Had this in here for multigene junctions, but going to do this in a later program;
run;


/* Cat together transcripts by junction */

proc sort data=xscript_junc_formatted;
    by junction_id transcript_id;
run;

/* get counts first */

proc freq noprint data=xscript_junc_formatted;
   tables junction_id / out=junc_count;
run;

proc sort data=junc_count;
  by descending count;
run;


*max=55 transcripts per junction;


data junctions_cat_xscript; 
  array xscripts[55] $ 53;

  retain xscripts1-xscripts55;

  set xscript_junc_formatted;
  by junction_id;
  
  if first.junction_id then do;
     call missing(of xscripts1-xscripts55);
     records = 0;
  end;

  records + 1;
  xscripts[records]=transcript_id;
  if last.junction_id then output;
run;

  *clean up the output file;

data junctions_cat_xscript2;
  set junctions_cat_xscript;
  length cat_xscript $ 2970;
  rename records= num_transcripts;
         cat_xscript= catx("|", OF xscripts1-xscripts55);
  drop xscripts1-xscripts55 transcript_id;
  rename cat_xscript=transcript_id;
  run;


/* Make permenant */

data splice.xscript_junctions;
   set junctions_cat_xscript2;
run;

/* remove temp datasets */

proc datasets noprint;
  delete TRANSCRIPT_JUNCTIONS xscript_junc_formatted junc_count
         junctions_cat_xscript junctions_cat_xscript2;
run;
quit;


/******** IMPORT SKIPPED EXON ANNOTATIONS  **********/


libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

/* Import skipped exon annotations */

    data WORK.EXON_SKIPPING_ANNOT    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_exon_skipping_annot.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat junction_id $79. ;
       informat flag_exonskip best32. ;
       informat num_skipped_exons best32. ;
       informat cat_skipped_exons $7489. ;
       format junction_id $79. ;
       format flag_exonskip best12. ;
       format num_skipped_exons best12. ;
       format cat_skipped_exons $7489. ;
    input
                junction_id $
                flag_exonskip
                num_skipped_exons
                cat_skipped_exons $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;


/*     data WORK.SKIPPED_EXON_LIST    ; */
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
/*     infile
 '/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_skipped_exons_list.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat junction_id $79. ;
        informat skipped_exon_id $39. ;
        informat flag_exonskip best32. ;
        format junction_id $79. ;
        format skipped_exon_id $39. ;
        format flag_exonskip best12. ;
     input
                 junction_id $
                 skipped_exon_id $
                 flag_exonskip
     ;
     if _ERROR_ then call symputx('_EFIERR_',1); */ /* set ERROR detection macro variable */
/*     run;  */


/* Make permenant */

data splice.EXON_SKIPPING_ANNOT;
  set EXON_SKIPPING_ANNOT;
  drop cat_skipped_exons; *don't need, but can reimport if necessary;
run;

/* if I need this later, I can reimport it
data splice.skipped_exon_list;
  set skipped_exon_list;
run; */

proc datasets noprint;
  delete EXON_SKIPPING_ANNOT skipped_exon_list;
run;
quit;


libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

/* Import IR events */


     data WORK.INTRON_RETENTION_EVENTS    ;
     %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
     infile
 '/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_intron_retention_176bp.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
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
 '/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_aceview_intron_retention_176bp.bed' delimiter='09'x MISSOVER DSD lrecl=32767 ;
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


/**** Merge logical junctions and annotated junctions */


libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

data xscript_junctions;
   set splice.xscript_junctions;
   rename junction_id=event_id;
   drop gene_id;
run;


/* Sort then merge */

proc sort data=splice.logical_junctions;
   by event_id;
run;

proc sort data=xscript_junctions;
   by event_id;
run;

/* Moment of truth - will they merge? */


data logical_junctions_w_xscript no_xscript no_logical;
    merge splice.logical_junctions (in=in1) xscript_junctions (in=in2);
    by event_id;
    if in1 and in2 then output logical_junctions_w_xscript;
    else if in1 then output no_xscript;
    else output no_logical;
run;


*12655954 logical junctions;
*608705 junctions from transcripts;
*608705 logical juncs with a xscript-junc match, yay!;
*12047249 logical juncs without a xscript-junc match;
*0 xscript juncs without a logical match, yay!!;

/* Okay so we know this works, need to output logical junctions with xscripts (if applicable) and a flag if annotated junction */

data logical_junctions_w_xscript no_logical_oops;
    merge splice.logical_junctions (in=in1) xscript_junctions (in=in2);
    by event_id;
    if in1 and in2 then do;
        flag_junction_annotated=1;
        output logical_junctions_w_xscript;
        end;
    else if in1 then do;
        num_transcripts=0;
        transcript_id=' ';
        flag_junction_annotated=0;
        output logical_junctions_w_xscript;
        end;
    else output no_logical_oops;
run;


/* Make permenant */

data splice.logical_junctions_w_xscript;
   set logical_junctions_w_xscript;
run;

/* Delete temp datasets */

proc datasets noprint;
  delete xscript_junctions no_xscript no_logical no_logical_oops logical_junctions_w_xscript;
run;
quit;



/***** Merge in exon skipping annotations *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

/* Sort and merge */

data exon_skipping_annot;
  set splice.exon_skipping_annot;
  rename junction_id=event_id;
run;


proc sort data=exon_skipping_annot;
   by event_id;
run;

proc sort data=splice.logical_junctions_w_xscript;
   by event_id;
run;


data junctions_w_exonskip no_exonskip no_junc_oops;
   merge splice.logical_junctions_w_xscript (in=in1) exon_skipping_annot (in=in2);
   by event_id;
   if in1 and in2 then output junctions_w_exonskip; *12655954 in, 12655954 out!;
   else if in1 then output no_exonskip; *0 obs, yay!;
   else output no_junc_oops; *0 obs, yay!;
run;

/* Make permenant */

data splice.junctions_w_exonskip;
   set junctions_w_exonskip;
run;

/* Delete temp datasets */

proc datasets noprint;
  delete junctions_w_exonskip exon_skipping_annot no_exonskip no_junc_oops;
run;
quit;


/***** Concatenate intron retention events *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';


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


/***** Add exon information *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';


/* Get donor/acceptor information to add to splicing events */
/* We want exon group, gene, flag_short_exon for donors and acceptors */

/* Get info from donors */

data donor_exon_info;
   set splice.donor_exons;
   drop donor_exon_start donor_exon_stop;
run;

/* Get info from acceptors */

data acceptor_exon_info;
   set splice.acceptor_exons;
   drop acceptor_exon_start acceptor_exon_stop;
run;


/* Sort donor and acceptor exon info */

proc sort data=donor_exon_info;
   by donor_exon;
run;

proc sort data=acceptor_exon_info;
   by acceptor_exon;
run;

/* Sort AS events by donor first */

proc sort data=splice.junction_and_ir_events;
   by donor_exon;
run;


/* Merge donor exon info - for now no_donor to check that all events without a donor exon are in fact IR events! */

data splicing_events_w_donors no_donor no_event;
   merge splice.junction_and_ir_events (in=in1) donor_exon_info (in=in2);
   by donor_exon;
   if in1 and in2 then output splicing_events_W_donors;
   else if in1 then output no_donor; * remove later if non-donors are ONLY introns!;
   else output no_event; *This can be more than zero, as not all exons will be used as donors! (ie, last exon per gene, single-exon genes);
run;

/* Check that all in no_donor are introns only */

proc freq data=no_donor noprint;
   tables donor_exon /out=no_donor_check;
run;

* all are introns!;


/* Sort AS events by acceptor first */

proc sort data=splice.junction_and_ir_events;
   by acceptor_exon;
run;


/* Merge donor exon info - for now no_donor to check that all events without a donor exon are in fact IR events! */

data splicing_events_w_acceptor no_acceptor no_event;
   merge splice.junction_and_ir_events (in=in1) acceptor_exon_info (in=in2);
   by acceptor_exon;
   if in1 and in2 then output splicing_events_W_acceptor;
   else if in1 then output no_acceptor; * remove later if non-acceptor are ONLY introns!;
   else output no_event; *This can be more than zero, as not all exons will be used as acceptor! (ie, last exon per gene, single-exon genes);
run;

/* Check that all in no_donor are introns only */

proc freq data=no_acceptor noprint;
   tables acceptor_exon /out=no_acceptor_check;
run;

* all are introns!;


/* Merge donor exon info - for now no_donor to check that all events without a donor exon are in fact IR events! */


proc sort data=splice.junction_and_ir_events;
   by donor_exon;
run;


data splicing_events_w_donors no_event;
   merge splice.junction_and_ir_events (in=in1) donor_exon_info (in=in2);
   by donor_exon;
   if in1 and in2 then output splicing_events_w_donors;
   else if in1 then do;
       donor_group=.;
       donor_gene='';
       flag_short_donor=.;
       flag_alt_donor=.;
       output splicing_events_w_donors;
       end;
   else output no_event;
run;


proc sort data=splicing_events_w_donors;
   by acceptor_exon;
run;


data splicing_events_w_acceptors no_event;
   merge splicing_events_w_donors (in=in1) acceptor_exon_info (in=in2);
   by acceptor_exon;
   if in1 and in2 then output splicing_events_w_acceptors;
   else if in1 then do;
       acceptor_group=.;
       acceptor_gene='';
       flag_short_acceptor=.;
       flag_alt_acceptor=.;
       output splicing_events_w_acceptors;
       end;
   else output no_event;
run;

/* Check genes */

data as_event_gene_check;
   set splicing_events_w_acceptors;
   if donor_gene = acceptor_gene then flag_bad_gene=0;
   else do;
       if donor_gene='' and acceptor_gene ne '' then flag_bad_gene=0;
       else if donor_gene ne '' and acceptor_gene = '' then flag_bad_gene=0;
       else flag_bad_gene=1;
       end;
run;


proc freq data=as_event_gene_check;
   tables flag_bad_gene;
run;
* all good!, collapse donor_gene and acceptor_gene into one variable: gene_id;
     

/* Make permenant */

data splice.splicing_events_w_exon_info;
   set as_event_gene_check;
   length gene_id $36.;
   if donor_exon='intron' then gene_id=acceptor_gene;
   else if acceptor_exon='intron' then gene_id=donor_gene;
   else gene_id=donor_gene;
   if flag_alt_donor=. then flag_alt_donor=0;
   if flag_alt_acceptor=. then flag_alt_acceptor=0;
   if flag_short_donor=. then flag_short_donor=0;
   if flag_short_acceptor=. then flag_short_acceptor=0;
   drop acceptor_gene donor_gene flag_bad_gene;
run;

proc datasets noprint;
  delete donor_exon_info acceptor_exon_info splicing_events_w_donors as_event_gene_check
   no_donor no_event splicing_events_w_acceptor splicing_events_w_acceptors no_acceptor;
run;
quit;




/***** Collapse duplicate events *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';
/* first count the number of "tiny" events - these are events that are smaller than the read size */
/* For the diabetes data this is 50bp */

data splicing_event_size;
   set splice.splicing_events_w_exon_info;
   event_size=donor_size+acceptor_size;
   if event_size lt 50 then flag_tiny_event=1;
   else flag_tiny_event=0;
   keep event_id event_size flag_tiny_event donor_size acceptor_size;
run;

proc sort data=splicing_event_size;
   by event_size donor_size acceptor_size;
run;

proc freq data=splicing_event_size noprint;
    tables flag_tiny_event / out=tiny_event_count;
run;

/* 9757360 events ge 150bp, 3176641 events lt 150bp */

 data small_events;
   set splicing_event_size;
   if event_size lt 150 then output;
run;

proc export data=small_events
   outfile='/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/small_exons_175bp.csv'
   dbms=csv
   replace;
run;

data splicing_events;
  set splice.splicing_events_w_exon_info;
  event_coords=cats(gene_id,"|",chr,":",donor_stop,":",acceptor_start,":",strand);
  event_size=donor_size+acceptor_size;
  keep event_id donor_exon acceptor_exon event_coords transcript_id  event_size num_transcripts;
run;

proc sort data=splicing_events;
  by event_coords event_size event_id;
run;

data first_event;
   set splicing_events;
   by event_coords;
   if first.event_coords then output;
   drop donor_exon acceptor_exon;
   rename event_id = new_event_id;
run;

/* use new_event_id to cat event_id, donor_exon, acceptor_exon */
proc sort data=first_event;
  by event_coords;
run;

proc sort data=splicing_events;
  by event_coords;
run;

data splicing_events_W_new_id oops1 oops2;
   merge first_event (in=in1) splicing_events (in=in2);
   by event_coords;
   if in1 and in2 then output splicing_events_w_new_id;
   else if in1 then output oops1;
   else output oops2;
run;

/* Add short exon info back in */
proc sort data=splicing_events_W_new_id;
  by new_event_id;
run;

proc sort data=splicing_event_size(rename=(event_id=new_event_id));
   by new_event_id;
run;

data collapsed_events_w_size event_not_ref no_size_oops;
   merge splicing_event_size (in=in1) splicing_events_W_new_id (in=in2);
   by new_event_id;
   if in1 and in2 then output collapsed_events_w_size;
   else if in1 then output event_not_ref;
   else output no_size_oops;
run;

data collapsed_events_for_cat;
   set collapsed_events_w_size;
   length donor_acceptor_num $13.;
   length donor_num $6.;
   length acceptor_num $6.;
   length gene_id $36.;
   if donor_exon='intron' then donor_num='intron';
   else donor_num=scan(donor_exon,2,':');
   if acceptor_exon='intron' then acceptor_num='intron';
   else acceptor_num=scan(acceptor_exon,2,':');
   donor_acceptor_num=catx('_', donor_num, acceptor_num);
   gene_id=scan(donor_exon,1,':');
   keep new_event_id donor_acceptor_num gene_id transcript_id donor_num acceptor_num;
run;

/* Clean up here */
proc datasets noprint;
delete splicing_event_size tiny_event_count small_events splicing_events first_event splicing_events_W_new_id oops1 oops2
collapsed_events_w_size event_not_ref no_size_oops;
run; quit;

/* Cat event_ids, donor_exons, acceptor_exons, transcript_ids; sum num_transcripts */
/*** cat event_ids ***/
data uniq_events;
   set collapsed_events_for_cat;
   keep new_event_id donor_acceptor_num;
   run;

proc sort data=uniq_events nodup;
   by new_event_id donor_acceptor_num;
run;

/* get counts first */
proc freq noprint data=uniq_events;
   tables new_event_id / out=event_count;
run;

proc sort data=event_count;
  by descending count;
run;
*max=377 events per ref_event;

data cat_event; 
  array event[377] $ 13.;
  retain event1-event377;
  set uniq_events;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of event1-event377);
     records = 0;
  end;
  records + 1;
  event[records]=donor_acceptor_num;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_event2;
  set cat_event;
  length cat_event_id $ 2500.;
  rename records= num_events;
         cat_event_id= catx("|", OF event1-event377);
  drop event1-event377 donor_acceptor_num;
  run;

/* Cat together donor exons */
data uniq_donor_exons;
   set collapsed_events_for_cat;
   keep new_event_id donor_num;
   run;

/* Cat together donors by event */
proc sort data=uniq_donor_exons nodup; *drop duplicated donor exons!;
    by new_event_id donor_num;
run;

/* get counts first */
proc freq noprint data=uniq_donor_exons;
   tables new_event_id / out=donor_count;
run;

proc sort data=donor_count;
  by descending count;
run;
*max=26 donor exons per event;

data cat_donors; 
  array donors[26] $ 6;
  retain donors1-donors26;
  set uniq_donor_exons;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of donors1-donors26);
     records = 0;
  end;
  records + 1;
  donors[records]=donor_num;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_donors2;
  set cat_donors;
  length donor_exons $ 182;
  rename records= num_donor_exons;
         donor_exons= catx("|", OF donors1-donors26);
  drop donors1-donors26 donor_num;
  run;

/* Cat together acceptor exons */
data uniq_acceptor_exons;
   set collapsed_events_for_cat;
   keep new_event_id acceptor_num;
   run;

/* Cat together acceptors by event */
proc sort data=uniq_acceptor_exons nodup; *drop duplicated acceptor exons!;
    by new_event_id acceptor_num;
run;

/* get counts first */
proc freq noprint data=uniq_acceptor_exons;
   tables new_event_id / out=acceptor_count;
run;

proc sort data=acceptor_count;
  by descending count;
run;
*max=29 acceptor exons per event;


data cat_acceptors; 
  array acceptors[29] $ 6;
  retain acceptors1-acceptors29;
  set uniq_acceptor_exons;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of acceptors1-acceptors29);
     records = 0;
  end;
  records + 1;
  acceptors[records]=acceptor_num;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_acceptors2;
  set cat_acceptors;
  length acceptor_exons $ 200;
  rename records= num_acceptor_exons;
         acceptor_exons= catx("|", OF acceptors1-acceptors29);
  drop acceptors1-acceptors29 acceptor_num;
  run;

/* Cat together transcripts */
data uniq_xscripts;
   set collapsed_events_for_cat;
   keep new_event_id transcript_id;
   run;

/* Cat together acceptors by event */
proc sort data=uniq_xscripts nodup; *drop duplicated transcripts!;
    by new_event_id transcript_id;
run;

/* get counts first */
proc freq noprint data=uniq_xscripts;
   tables new_event_id / out=xscripts_count;
run;

proc sort data=xscripts_count;
  by descending count;
run;
*max=29 transcripts per event (actually more, but this is 29 "concatenated transcripts" per event);

/* Commenting out for now. This is making a 200GB+ SAS dataset! */

 data cat_xscripts; 
  array xscripts[29] $ 1085;
  retain xscripts1-xscripts29;
  set uniq_xscripts;
  by new_event_id;
  if first.new_event_id then do;
     call missing(of xscripts1-xscripts29);
     records = 0;
  end;
  records + 1;
  xscripts[records]=transcript_id;
  if last.new_event_id then output;
run;

  *clean up the output file;
data cat_xscripts2;
  set cat_xscripts;
  length transcripts $ 2000;
         transcripts= catx("|", OF xscripts1-xscripts29);
  drop xscripts1-xscripts29 transcript_id records;
  rename transcripts=transcript_id;
  run;
*/
/* Count number of unique transcripts per collapsed event */

data uniq_xscript_counts;
   set uniq_xscripts;
   if transcript_id='' then xscript_count=0;
   else xscript_count=count(transcript_id,'|')+1;
   keep new_event_id xscript_count;
run;


/* Sum counts */
proc sort data=uniq_xscript_counts;
   by new_event_id;
run;

proc means data=uniq_xscript_counts noprint;
   by new_event_id;
   var xscript_count;
   output out=num_xscripts sum=total_transcripts;
run;

/* Merge together */
proc sort data=cat_event2;
by new_event_id;
run;

/*proc sort data=cat_xscripts2;
by new_event_id;
run;*/

proc sort data=cat_donors2;
by new_event_id;
run;

proc sort data=cat_acceptors2;
by new_event_id;
run;

proc sort data=num_xscripts;
by new_event_id;
run;

data collapsed_event_info;
   *merge cat_event2 num_xscripts cat_xscripts2 cat_donors2 cat_acceptors2;
   merge cat_event2 num_xscripts cat_donors2 cat_acceptors2;
   by new_event_id;
   drop _TYPE_ _FREQ_;
   rename new_event_id=event_id;
run;

/* Collapse full event information */
data splice_events_all;
   set splice.splicing_events_w_exon_info;
   drop donor_exon acceptor_exon num_transcripts transcript_id;
run;

proc sort data=collapsed_event_info;
   by event_id;
run;

proc sort data=splice_events_all;
   by event_id;
run;

data collapsed_events_final not_ref no_event_oops;
   merge collapsed_event_info (in=in1) splice_events_all (in=in2);
   by event_id;
   if in1 and in2 then output collapsed_events_final;
   else if in1 then output no_event_oops;
   else output not_ref;
run;

/* Make permenant */
data splice.collapsed_splicing_events;
    length new_event_id $2550.;
    length feature1_id $250.;
    length feature2_id $250.;
    set collapsed_events_final;
    new_event_id=catx(':', gene_id, cat_event_id);
    feature1_id=catx(':', gene_id, donor_exons);
    feature2_id=catx(':', gene_id, acceptor_exons);
    rename total_transcripts=num_transcripts;
    event_size=donor_size+acceptor_size;
    if event_size lt 50 then flag_tiny_event=1;
    else flag_tiny_event=0;
    drop event_id cat_event_id;
run;

/* Delete temp datasets */
proc datasets noprint;
delete splicing_event_size tiny_event_count small_events
splicing_events first_event splicing_events_W_new_id oops1 oops2
collapsed_events_w_size event_not_ref no_size_oops
collapsed_events_for_cat uniq_events event_count
cat_event cat_event2 uniq_donor_exons
donor_count cat_donors cat_donors2 uniq_acceptor_exons
acceptor_count cat_acceptors cat_acceptors2 uniq_xscripts
xscripts_count cat_xscripts cat_xscripts2 uniq_xscript_counts
num_xscripts collapsed_event_info splice_events_all
collapsed_events_final not_ref no_event_oops;
run; quit;





/***** Format splicing event annotations *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';

data splicing_events_formatting;
   set splice.collapsed_splicing_events;
   length feature1_type $13.;
   length feature2_type $13.;
   rename new_event_id=event_id;
   if strand='+' then do;
       if event_type='exon_junction' then do;
            feature1_type='exon_donor';
            feature2_type='exon_acceptor';
            end;
       else if event_type='intron_retention' then do;
            feature1_type='exon_donor';
            feature2_type='intron';
            end;
       end;
   else if strand='-' then do; 
       if event_type='exon_junction' then do;
            feature1_type='exon_acceptor';
            feature2_type='exon_donor';
            end;
       else if event_type='intron_retention' then do;
            feature1_type='intron';
            feature2_type='exon_donor';
            end;
       end;
run;

/* Rearrange columns */


data splicing_events_reordered;
   retain event_id event_type num_events gene_id chr strand event_size flag_tiny_event num_transcripts transcript_id feature1_id feature1_type
feature2_id feature2_type donor_start donor_stop num_donor_exons donor_exons donor_size flag_short_donor acceptor_start acceptor_stop num_acceptor_exons
acceptor_exons acceptor_size flag_short_acceptor flag_junction_annotated flag_intron_retention num_skipped_exons flag_exonskip flag_alt_donor
flag_alt_acceptor;
   set splicing_events_formatting;
rename flag_tiny_event=flag_event_short donor_start=feature1_start donor_stop=feature1_stop num_donor_exons=num_feature1 donor_exons=feature1_list donor_size=feature1_size flag_short_donor=flag_feature1_short acceptor_start=feature2_start acceptor_stop=feature2_stop num_acceptor_exons=num_feature2 acceptor_exons=feature2_list acceptor_size=feature2_size flag_short_acceptor=flag_feature2_short;
drop acceptor_group donor_group;
run;


/* Make permenant */

data splice.splicing_events_annotations;
    set splicing_events_reordered;
run;


/* Delete temp datasets */

proc datasets noprint;
  delete splicing_events_formatting splicing_events_reordered;
run;quit;

/* Export catalogue CSV */

proc export data=splice.splicing_events_annotations
    outfile='/home/jrbnewman/McLab/junction_annotations/pipeline_output/aceview_hg19/hg19_splicing_catalogue_176bp.csv'
    dbms=csv
    replace;
run;

58860600
58860564
      36

/***** Catalogue2BED *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';


/* Update event lengths to be 74bp */

data splice_new_lengths;
  set splice.splicing_events_annotations;

  /* Update features */
  if feature1_type='intron' then do;
       feature1_size2=63;
       feature1_start2=feature1_stop-63;
       end;
  else if feature1_size lt 63 then do;
       feature1_size2=feature1_size;
       feature1_start2=feature1_start;
       end;
  else do;
      feature1_size2=63;
      feature1_start2=feature1_stop-63;
      end;

  if feature2_type='intron' then do;
      feature2_size2=63;
      feature2_stop2=feature2_start+63;
      end;
  else if feature2_size lt 63 then do;
      feature2_size2=feature2_size;
      feature2_stop2=feature2_stop;
      end;
  else do;
      feature2_size2=63;
      feature2_stop2=feature2_start+63;
      end;
  event_size2=feature1_size2+feature2_size2;
run;





data info_for_bed;
   length score $1.;
   length color $7.;
   set splice_new_lengths;
   score='.';
   color='255,0,0';
   block1_start=0;

   if event_type='intron_retention' then do;
      num_blocks=1;
      if strand='+' then do;
         block1_length=event_size2+1;
         totalstart1=feature1_start2;
         totalstop1=feature2_stop2+1;
         totalstart2=feature1_start2;
         totalstop2=feature2_stop2+1;
         end;
      if strand='-' then do;
         block1_length=event_size2+1;
         totalstart1=feature1_start2-1;
         totalstop1=feature2_stop2;
         totalstart2=feature1_start2-1;
         totalstop2=feature2_stop2;
         end;
      end;
   else do;
      num_blocks=2;
          totalstart1=feature1_start2;
          totalstop1=feature2_stop2;
          totalstart2=feature1_start2;
          totalstop2=feature2_stop2;
      block1_length=feature1_size2;
      block2_length=feature2_size2;
      block2_start=feature2_start-feature1_start2;
      end;
   keep event_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand score color num_blocks block1_length block2_length block1_start block2_start;
run;


/* BED12 format 


chrom		chr
totalStart	feature1_start
totalStop	feature2_stop
name		event_id
score		"."
strand		strand
totalStart	feature1_start
totalStop	feature2_stop
color		255,0,0
num		"1" if intron retention, "2" if junction
lengths		block lengths (total if intron_retention, feature1_stop-feature1_start,feature2_stop-feature2_start
starts		if retention then 0, if junction then (0,feature2_start-feature1_start)

*/

data assemble_bed;
    retain chr;
    retain totalstart1;
    retain totalstop1;
    retain event_id;
    retain score;
    retain strand;
    retain totalstart2;
    retain totalstop2;
    retain color;
    retain num_blocks;
    length lengths_cat $10.;
    length starts_cat $10.;
    set info_for_bed;
    if num_blocks=2 then do;
        lengths_cat=catx(',', block1_length, block2_length);
        starts_cat=catx(',', block1_start, block2_start);
        end;
    else if num_blocks=1 then do;
        lengths_cat=put(block1_length, 3.);
        starts_cat=put(block1_start, 1.);
        end;
    drop block1_length block2_length block1_start block2_start;
run;

proc sort data=assemble_bed;
   by chr totalstart1 totalstop1 strand;
run;


/* output as BED */

proc export data=assemble_bed
	outfile='/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/hg19_splicing_catalogue_126bp_redo.bed'
	dbms=tab replace;
        putnames=no;
	run;


proc datasets noprint;
   delete assemble_bed info_for_bed;
run;
quit;



/* Update to be 126 bp */


data splice_new_lengths;
  set splice.splicing_events_annotations;

  /* Update features */
  if feature1_size lt 63 then do;
      feature1_size2=feature1_size;
      feature1_start2=feature1_start;
      end;
  else do;
      feature1_size2=63;
      feature1_start2=feature1_stop-63;
      end;
  if feature2_size lt 63 then do;
      feature2_size2=feature2_size;
      feature2_stop2=feature2_stop;
      end;
  else do;
      feature2_size2=63;
      feature2_stop2=feature2_start+63;
      end;
  event_size2=feature1_size2+feature2_size2;
run;





data info_for_bed;
   length score $1.;
   length color $7.;
   set splice_new_lengths;
   score='.';
   color='255,0,0';
   block1_start=0;

   if event_type='intron_retention' then do;
      num_blocks=1;
      if strand='+' then do;
         block1_length=event_size2+1;
         totalstart1=feature1_start2;
         totalstop1=feature2_stop2+1;
         totalstart2=feature1_start2;
         totalstop2=feature2_stop2+1;
         end;
      if strand='-' then do;
         block1_length=event_size2+1;
         totalstart1=feature1_start2-1;
         totalstop1=feature2_stop2;
         totalstart2=feature1_start2-1;
         totalstop2=feature2_stop2;
         end;
      end;
   else do;
      num_blocks=2;
          totalstart1=feature1_start2;
          totalstop1=feature2_stop2;
          totalstart2=feature1_start2;
          totalstop2=feature2_stop2;
      block1_length=feature1_size2;
      block2_length=feature2_size2;
      block2_start=feature2_start-feature1_start2;
      end;
   keep event_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand score color num_blocks block1_length block2_length block1_start block2_start;
run;


/* BED12 format 


chrom		chr
totalStart	feature1_start
totalStop	feature2_stop
name		event_id
score		"."
strand		strand
totalStart	feature1_start
totalStop	feature2_stop
color		255,0,0
num		"1" if intron retention, "2" if junction
lengths		block lengths (total if intron_retention, feature1_stop-feature1_start,feature2_stop-feature2_start
starts		if retention then 0, if junction then (0,feature2_start-feature1_start)

*/

data assemble_bed;
    retain chr;
    retain totalstart1;
    retain totalstop1;
    retain event_id;
    retain score;
    retain strand;
    retain totalstart2;
    retain totalstop2;
    retain color;
    retain num_blocks;
    length lengths_cat $10.;
    length starts_cat $10.;
    set info_for_bed;
    if num_blocks=2 then do;
        lengths_cat=catx(',', block1_length, block2_length);
        starts_cat=catx(',', block1_start, block2_start);
        end;
    else if num_blocks=1 then do;
        lengths_cat=put(block1_length, 3.);
        starts_cat=put(block1_start, 1.);
        end;
    drop block1_length block2_length block1_start block2_start;
run;

proc sort data=assemble_bed;
   by chr totalstart1 totalstop1 strand;
run;


/* output as BED */

proc export data=assemble_bed
	outfile='/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/hg19_splicing_catalogue_126bp_redo.bed'
	dbms=tab replace;
        putnames=no;
	run;


proc datasets noprint;
   delete assemble_bed info_for_bed;
run;
quit;







