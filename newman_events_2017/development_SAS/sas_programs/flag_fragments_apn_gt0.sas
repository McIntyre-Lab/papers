libname events '!MCLAB/event_analysis/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* For each cell type, flag if fragments are detected */

data set_group;
  length cell_type $3.;
  set events.mm10_refseq_fragment_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type fragment_id apn;
run;

data flag_frag_apn_gt0;
  set set_group;
  if apn > 0 then flag_fragment_apn_gt0=1;
  else flag_fragment_apn_gt0=0;
run;


proc sort data=flag_frag_apn_gt0;
  by cell_type fragment_id;
proc means data=flag_frag_apn_gt0 noprint;
  by cell_type fragment_id;
  var flag_fragment_apn_gt0;
  output out=mean_on mean=mean_gt0;
run;

data nsc old;
  set mean_on;
  if cell_type="NSC" then output nsc;
  else output old;
run;

data flag_on_nsc;
  set nsc;
  if mean_gt0 ge 0.5 then flag_fragment_nsc_on=1;
  else flag_fragment_nsc_on=0;
  keep fragment_id flag_fragment_nsc_on;
run;

data flag_on_old;
  set old;
  if mean_gt0 ge 0.5 then flag_fragment_old_on=1;
  else flag_fragment_old_on=0;
  keep fragment_id flag_fragment_old_on;
run;

proc sort data=flag_on_nsc;
   by fragment_id;
proc sort data=flag_on_old;
   by fragment_id;
run;

data events.flag_fragment_on;
   merge flag_on_nsc (in=in1) flag_on_old (in=in2);
   by fragment_id;
   if in1 and in2;
run;

