libname events '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* For each cell type, flag if fragments are detected */

data set_group;
  length cell_type $3.;
  set event.mm10_refseq_intron_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type intron_id apn;
run;

data flag_intron_apn_gt0;
  set set_group;
  if apn > 0 then flag_intron_apn_gt0=1;
  else flag_intron_apn_gt0=0;
run;


proc sort data=flag_intron_apn_gt0;
  by cell_type intron_id;
proc means data=flag_intron_apn_gt0 noprint;
  by cell_type intron_id;
  var flag_intron_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data nsc old;
  set mean_on;
  if cell_type="NSC" then output nsc;
  else output old;
run;

data flag_on_nsc;
  set nsc;
  if mean_gt0 ge 0.5 then flag_intron_nsc_on=1;
  else flag_intron_nsc_on=0;
  keep intron_id flag_intron_nsc_on;
run;

data flag_on_old;
  set old;
  if mean_gt0 ge 0.5 then flag_intron_old_on=1;
  else flag_intron_old_on=0;
  keep intron_id flag_intron_old_on;
run;

proc sort data=flag_on_nsc;
   by intron_id;
proc sort data=flag_on_old;
   by intron_id;
run;

data event.flag_intron_on;
   merge flag_on_nsc (in=in1) flag_on_old (in=in2);
   by intron_id;
   if in1 and in2;
run;

