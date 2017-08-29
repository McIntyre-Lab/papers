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
