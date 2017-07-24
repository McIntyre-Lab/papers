ods listing; ods html close;

libname event "!MCLAB/event_analysis/sas_data";
libname mm10 "!MCLAB/useful_mouse_data/mm10/sas_data";


/* Count unannotated junctions by type and by whether they have a PB hit or not */

data junc_w_hit;
  set event.blast_junc_to_pb_w_annot;
  where flag_intron_retention=0 and flag_feature_on_ge5=1
        and flag_event_has_hit=1;
  keep event_id;
run;

data unannot_junc;
  set event.unannot_junc_dtct_flag_fus;
  keep event_id ;
run;


proc sort data=junc_w_hit nodup;
   by event_id;
proc sort data=unannot_junc;
   by event_id;
run;

data junc_pb_hit;
  merge unannot_junc (in=in1) junc_w_hit (in=in2);
  by event_id;
  if in2 then flag_pb_hit=1; else flag_pb_hit=0;
  if in1 then output;
run;


proc freq data=junc_pb_hit;
   tables flag_pb_hit;
run;

proc print data=ir_count;
run;


/*

                                         Cumulative    Cumulative
 flag_pb_hit    Frequency     Percent     Frequency      Percent
 ----------------------------------------------------------------
           0         182       31.22           182        31.22
           1         401       68.78           583       100.00
                                                 
*/


