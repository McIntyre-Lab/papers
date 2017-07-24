ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Check that the set of genes with DE MISO SE are the same set with DE/DD SE in Event analysis */

/* Run DE analysis -- basic model, just look at log_APN vs group */

* Get events for analysis;
  
data events_for_analysis;
   set event.flag_splicing_on;
   where flag_event_nsc_on=1 and flag_event_old_on=1;
   keep event_id;
run; *136131 events analyzable;

* get exonskipping events;
data exonskip;
   set evspl.splicing_events_annot_refseq;
   where flag_exonskip=1;
   keep event_id gene_id;
run;

proc sort data=exonskip;
  by event_id;
proc sort data=events_for_analysis;
   by event_id;
run;

data events_for_analysis2;
   merge events_for_analysis (in=in1) exonskip (in=in2);
   by event_id;
   if in1 and in2;
run; *13907 exonskip events analyzable;

data set_group;
  length cell_type $3.;
  set event.mm10_refseq_splicing_counts;
  if sample_id in ('NSC1','NSC2') then cell_type="NSC";
  if sample_id in ('OLD1','OLD2') then cell_type="OLD";
  keep sample_id cell_type event_id apn;
run;

proc sort data=set_group;
   by event_id;
run;

data counts_for_anova;
   merge events_for_analysis2 (in=in1) set_group (in=in2);
   by event_id;
   if in1 and in2;
   log_apn = log(apn + 1);
run;

proc sort data=counts_for_anova;
   by event_id cell_type;
run;



ods listing close;
proc glimmix data=counts_for_anova;
   by event_id;
   class cell_type;
   model log_apn = cell_type / htype=1;
   output out=resid resid=resid pred=pred student=stu;
ods output tests1=anova;
run;
quit;

/* Flag residuals */
proc univariate data=resid normal noprint;
   by event_id;
   var Resid;
   output out=normtest probn=pnorm;
run;

data flag_resids;
  set normtest;
  if pnorm = . then flag_fail_norm=.;
  else if pnorm le 0.05 then flag_fail_norm=1;
  else flag_fail_norm=0;
run;

proc freq data=flag_Resids noprint;
  tables flag_fail_norm / out=splicing_flag_fail_norm;
run;

ods listing;
proc print data=splicing_flag_fail_norm;
run;

/*
        flag_
        fail_
 Obs     norm    COUNT    PERCENT

  1       .          6      .
  2       0      12500    89.9216
  3       1       1401    10.0784
*/

/* Make permenant */

data event.exonskip_flag_failnorm;
   set splicing_flag_fail_norm;
run;

data event.exonskip_anova;
  set anova;
run;

data event.exonskip_flag_resid;
   set flag_resids;
run;

