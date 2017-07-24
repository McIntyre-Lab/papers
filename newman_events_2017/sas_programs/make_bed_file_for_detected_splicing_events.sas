
libname event '!MCLAB/event_analysis/sas_data';
libname splice '!MCLAB/conesa_pacbio/sas_data/splicing';

/* Create a BED file for the set of detected splicing events so that I can extract the FASTA sequences and BLAST onto
   the set of PacBio transcripts */

data events_on;
  set event.splicing_on_apn_gt0;
  where  flag_splicing_on=1;
  keep event_id;
run;

data event2coord;
   set splice.splicing_events_annot_refseq;
run;

proc sort data=event2coord;
   by event_id;
proc sort data=events_on;
  by event_id;
run;

data events_for_bed;
   merge events_on (in=in1) event2coord (in=in2);
   by event_id;
   if in1 and in2;
run;

data info_for_bed;
   length score $1.;
   length color $7.;
   set events_for_bed;
   score='.';
   color='255,0,0';
   block1_start=0;

   if event_type='intron_retention' then do;
      num_blocks=1;
      if strand='+' then do;
         block1_length=event_size+1;
         totalstart1=feature1_start;
         totalstop1=feature2_stop+1;
         totalstart2=feature1_start;
         totalstop2=feature2_stop+1;
         end;
      if strand='-' then do;
         block1_length=event_size+1;
         totalstart1=feature1_start-1;
         totalstop1=feature2_stop;
         totalstart2=feature1_start-1;
         totalstop2=feature2_stop;
         end;
      end;
   else do;
      num_blocks=2;
          totalstart1=feature1_start;
          totalstop1=feature2_stop;
          totalstart2=feature1_start;
          totalstop2=feature2_stop;
      block1_length=feature1_size;
      block2_length=feature2_size;
      block2_start=feature2_start-feature1_start;
      end;
   keep event_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand score color num_blocks block1_length block2_length block1_start block2_start;
run;


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
        outfile='!MCLAB/event_analysis/references/refseq_mm10_detected_splicing_events_nsc.bed'
        dbms=tab replace;
        putnames=no;
        run;

