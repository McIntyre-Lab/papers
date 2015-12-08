/***** Catalogue2BED *******/

libname splice '/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/splice';


data info_for_bed;
   length score $1.;
   length color $7.;
   set splice.splicing_events_annotations;
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
	outfile='/media/jrbnewman/ac89a883-cbf2-4ed0-8e2f-19c25fead575/hg19_splicing_catalogue_74bp.bed'
	dbms=tab replace;
        putnames=no;
	run;


proc datasets noprint;
   delete assemble_bed info_for_bed;
run;
quit;

