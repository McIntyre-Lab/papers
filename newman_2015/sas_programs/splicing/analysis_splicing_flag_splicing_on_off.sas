
/********************* FLAG COUNTS FOR SPLICING EVENTS *************************/

/* library */
libname mysas '/scratch/lfs/patcon/jnewman/sas_analysis/sas_data2';


/* Add in flags for counts etc. */



data design_file_&sysparm.;
    set mysas.design_file_se;
    keep Name subject_id cell_type;
run;


proc sort data=design_file_&sysparm. nodup;
  by Name;
  run;

proc sort data=mysas.mapped_reads_summed_byevent_&sysparm.;
      by Name;
run;

 data mapped_reads_summed_&sysparm.;
   merge mysas.mapped_reads_summed_byevent_&sysparm. (in=in1) design_file_&sysparm. (in=in2);
   by Name;
   *logapn = log(apn +1);
   if in1;
 run;


data counts_flag_on_&sysparm.;
   set mapped_reads_summed_&sysparm.;
   if depth ge 10 then flag_event_on=1;
   if flag_event_on=1 then output;
   keep event_id cell_type flag_event_on;
run;


proc sort data=counts_flag_on_&sysparm.;
   by event_id cell_type flag_event_on;
run;

proc freq data=counts_flag_on_&sysparm. noprint;
   by event_id cell_type;
   tables flag_event_on / out=counts_subjects_exp_&sysparm.;
run;

*want a flag_CD19_on flag_CD4_on flag_CD8_on flag_on_all;

data counts_on_only_&sysparm.;
   set counts_subjects_exp_&sysparm.;
   if flag_event_on=1 then output;
   drop percent;
run;

data mysas.counts_on_only_&sysparm.;
   set counts_on_only_&sysparm.;
run;



proc sort data=counts_on_only_&sysparm.;
  by event_id;
run;

data CD19_on_&sysparm. CD4_on_&sysparm. CD8_on_&sysparm.;
   set counts_on_only_&sysparm.;
   if cell_type = 'CD19' then output CD19_on_&sysparm.;
   if cell_type = 'CD4' then output CD4_on_&sysparm.;
   if cell_type = 'CD8' then output CD8_on_&sysparm.;
run;

data CD19_on_2_&sysparm.;
    set CD19_on_&sysparm.;
    if count>40 then flag_CD19_on=1;
    else flag_CD19_on=0;
    rename count=count_CD19;
    drop cell_type;
run;

data CD4_on_2_&sysparm.;
    set CD4_on_&sysparm.;
    if count>40 then flag_CD4_on=1;
    else flag_CD4_on=0;
    rename count=count_CD4;
    drop cell_type;
run;


data CD8_on_2_&sysparm.;
    set CD8_on_&sysparm.;
    if count>40 then flag_CD8_on=1;
    else flag_CD8_on=0;
    rename count=count_CD8;
    drop cell_type;
run;

proc sort data=CD19_on_2_&sysparm.;
   by event_id;
run;

proc sort data=CD4_on_2_&sysparm.;
   by event_id;
run;

proc sort data=CD8_on_2_&sysparm.;
   by event_id;
run;



proc sort data=mapped_reads_summed_&sysparm.;
   by event_id;
run;

data flags_on_&sysparm.;
   merge CD19_on_2_&sysparm. CD4_on_2_&sysparm. CD8_on_2_&sysparm.;
   by event_id;
   if flag_CD19_on=. then flag_CD19_on=0;
   if flag_CD4_on=. then flag_CD4_on=0;
   if flag_CD8_on=. then flag_CD8_on=0;
   if flag_CD19_on=1 and flag_CD4_on=1 and flag_CD8_on=1 then flag_all_on=1;
   else flag_all_on=0;
run;


data mysas.counts_by_splicing_w_flags_&sysparm. oops no_counts;
   merge mapped_reads_summed_&sysparm. (in=in1) flags_on_&sysparm. (in=in2);
   by event_id;
   if in1 and in2 then output mysas.counts_by_splicing_w_flags_&sysparm.;
   else if in1 then output oops;
   else output no_counts;
run;

