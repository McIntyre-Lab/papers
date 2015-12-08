/* Flagging logical junctions if they are annotated to a transcript */
/* libraries */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';

/* import all possible, logical exon junctions */


/*proc import datafile='/home/jrbnewman/McLab/junction_annotations/generated_files/hg19_logical_junctions_info.csv'
    out=junctions_logical
    dbms=csv
    replace;
    getnames=yes;
    guessingrows=2700000;
run;*/


    data WORK.JUNCTIONS_LOGICAL    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
'/home/jrbnewman/McLab/junction_annotations/generated_files/hg19_logical_junctions_info.csv'
delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat chr $2. ;
       informat donor_start best32. ;
       informat donor_stop best32. ;
       informat acceptor_start best32. ;
       informat acceptor_stop best32. ;
       informat junction_id $86. ;
       informat strand $1. ;
       informat exonA $61. ;
       informat exonB $61. ;
       format chr $2. ;
       format donor_start best12. ;
       format donor_stop best12. ;
       format acceptor_start best12. ;
       format acceptor_stop best12. ;
       format junction_id $86. ;
       format strand $1. ;
       format exonA $61. ;
       format exonB $61. ;
    input
                chr $
                donor_start
                donor_stop
                acceptor_start
                acceptor_stop
                junction_id $
                strand $
                exonA $
                exonB $6400 Southwest 20th Ave
    ;
 if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
 run;

/* Sort */


proc sort data=splice.junctions_xscript_annotated_hg19 out=junctions_w_xscript;
   by junction_id;
run;

proc sort data=junctions_logical;
   by chr donor_stop acceptor_start strand;
run;

/* Add in variable to match by - chr:start:stop:strand */

data junctions_xscripts_tag;
   set junctions_w_xscript;
   rename junction_id=junction_pos;
run;

data junctions_logical_tag;
   set junctions_logical;
   length junction_pos $50.;
   junction_pos=catx(':',chr,donor_stop,acceptor_start,strand);
run;

/* sort by junction_pos */

proc sort data=junctions_logical_tag;
   by junction_pos;
run;

proc sort data=junctions_xscripts_tag;
   by junction_pos;
run;


/* merge and flag_annotated */
data junc_flag_annotated no_logical oops;
   merge junctions_logical_tag (in=in1) junctions_xscripts_tag (in=in2);
   by junction_pos;
   if in1 and in2 then do;
      flag_junction_annotated=1;
      output junc_flag_annotated;
      end;
   else if in1 then do;
      flag_junction_annotated=0;
      output junc_flag_annotated;
      end;
   else if in2 then output no_logical;
   else output oops; *0 obs, woo!;
run;

*2610258 logical junctions, 199181 annotated junctions;
*64062 junctions with no logical counterpart;
*probably due to the way logical junctions are generated;
*and we have genes and xscripts for each annotated junction;


/* quick count of annotated vs non-annotated junctions */

proc freq data=junc_flag_annotated;
   tables flag_junction_annotated;
run;

*2346890 junctions are not annotated to a transcript;
*263368 junctions annotated to a transcript;

/* Make datasets permenant */

data splice.junctions_annotated_hg19;
    set junc_flag_annotated;
run;


/* Clean up */
    proc datasets nolist;
        delete junctions_logical junctions_xscripts_tag junctions_logical_tag
        junctions_w_xscript junc_flag_annotated no_logical oops
        ;
        run;
        quit;


