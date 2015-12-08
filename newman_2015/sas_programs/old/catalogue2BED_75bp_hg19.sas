/* Getting BED info */
/* libraries */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';
libname splice2 '/media/jrbnewman/SAS_WRK1/';

data info_for_bed;
   length score $1.;
   length color $7.;
   set splice2.splicing_event_catalogue_hg19;
   score='.';
   color='255,0,0';
   block1_start=0;

   if flag_intron_retention=1 then do;
      num_blocks=1;
      if strand="+" then do;
          block1_length=feature2_stop-feature1_start;
          totalstop2=feature2_stop;
          totalstart2=feature1_start;
          totalstop1=feature2_stop;
          totalstart1=feature1_start;
          other_stop=feature2_stop;
          other_start=feature1_start;
          end;
      else do;
          block1_length=feature1_start-feature2_stop;
          totalstart1=feature2_stop;
          totalstop1=feature1_start;
          totalstart2=feature2_stop;
          totalstop2=feature1_start;
          other_start=feature2_stop;
          other_stop=feature1_start;
          end;
      end;
   else do;
      num_blocks=2;
          totalstop2=feature2_stop;
          totalstart2=feature1_start;
          totalstop1=feature2_stop;
          totalstart1=feature1_start;
      block1_length=feature1_stop-feature1_start;
      block2_length=feature2_stop-feature2_start;
      block2_start=feature2_start-feature1_start;
      end;
   keep event_id chr totalstart1 totalstop1 totalstart2 totalstop2 strand flag_intron_retention score color num_blocks block1_length block2_length block1_start block2_start other_start other_stop;
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
    retain totalstop2;
    retain event_id;
    retain score;
    retain strand;
    retain other_start;
    retain other_stop;
    retain color;
    retain num_blocks;
    length lengths_cat $10.;
    length starts_cat $10.;
    set info_for_bed;
    other_start=totalstart1;
    other_stop=totalstop2;
    if num_blocks=2 then do;
        lengths_cat=catx(',', block1_length, block2_length);
        starts_cat=catx(',', block1_start, block2_start);
        end;
    else if num_blocks=1 then do;
        lengths_cat=put(block1_length, 3.);
        starts_cat=put(block1_start, 1.);
        end;
    drop flag_intron_retention block1_length block2_length block1_start block2_start totalstop1 totalstart2;
run;



/* output as BED */

proc export data=assemble_bed
	outfile='/home/jrbnewman/McLab/junction_annotations/pipeline_output/hg19_splicing_catalogue_75bp.bed'
	dbms=tab replace;
        putnames=no;
	run;

