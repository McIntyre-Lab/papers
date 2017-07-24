ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* For detected unannotated junction, flag if the adjacent 5' and 3' fragments and fusions are detected at APN>5 (both reps) */

data donor_fus_on;
   set event.flag_fusion_on_apn5;
   keep fusion_id flag_fusion_nsc_on;
   rename fusion_id=donor_fusion_id flag_fusion_nsc_on=flag_donor_on_ge5;
run;

data acceptor_fus_on;
   set event.flag_fusion_on_apn5;
   keep fusion_id flag_fusion_nsc_on;
   rename fusion_id=acceptor_fusion_id flag_fusion_nsc_on=flag_acceptor_on_ge5;
run;

data event_on;
   set event.flag_splicing_on_apn5;
   keep event_id flag_event_nsc_on;
   rename flag_event_nsc_on=flag_event_on_ge5;
run;

data unannot_junc2fus;
   set event.unannot_jnc_w_flanking_fus;
run;

proc sort data=unannot_junc2fus;
   by event_id;
proc sort data=event_on;
   by event_id;
run;

data unannot_flag_event;
  merge unannot_junc2fus (in=in1) event_on (in=in2);
  by event_id;
  if in1 and in2;
run;

proc sort data=unannot_flag_event;
   by donor_fusion_id;
proc sort data=donor_fus_on;
   by donor_fusion_id;
run;

data unannot_flag_donor;
   merge unannot_flag_event (in=in1) donor_fus_on (in=in2);
   by donor_fusion_id;
   if in1 and in2;
run;

proc sort data=unannot_flag_donor;
   by acceptor_fusion_id;
proc sort data=acceptor_fus_on;
   by acceptor_fusion_id;
run;

data unannot_flag_acceptor;
   merge unannot_flag_donor (in=in1) acceptor_fus_on (in=in2);
   by acceptor_fusion_id;
   if in1 and in2;
run;

/* Crosstabs on donor*acceptor*junction on */

proc freq data=unannot_flag_acceptor noprint;
  where flag_event_on_ge5=1;
  tables flag_donor_on_ge5*flag_acceptor_on_ge5 / out=unannot_flag_check;
run;

proc print data=unannot_flag_check;
run;

/*


  flag_
 donor_    flag_acceptor_
 on_ge5        on_ge5        COUNT

    .             0             1
    .             1            11
    0             .             1
    0             0             9
    0             1            17
    1             .            11
    1             0            10
    1             1           583

583 junctions to keep!
*/

/* Subset detected unannotated junctions and make permenant */

data event.unannot_junc_dtct_flag_fus;
   set unannot_flag_acceptor;
   if flag_event_on_ge5=1 and flag_donor_on_ge5=1 and flag_acceptor_on_ge5=1;
run;

