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

