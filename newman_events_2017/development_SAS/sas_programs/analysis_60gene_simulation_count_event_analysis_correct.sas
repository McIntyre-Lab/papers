/* Libraries */
libname event  '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* 60 gene simulation, Event Analysis check:
   what do we get right in terms of transcript retention?

   I will look at: 100% APN0, 75% APN0, 75% APN5 */

/* Get list of transcripts simulated */

data sim_list;
  set event.polyester_xs_list_60genes;
run;

/* For each of the 3 reduced references count the overlap between the simulated transcripts
   and the ones selected by EA:
   (1) Correct transcripts
   (2) Additional, unrelated transcripts (no "related" since we simulate all transcripts from
       these genes
   (3) Missing transcripts  */

%macro countEA(propDtct,apnLvl);

data ea_list;
   set event.bin_xs_by_dtct_apn&apnLvl._10gn;
   where perc_features_dtct ge &propDtct.;
run;

proc sort data=ea_list;
  by transcript_id;
proc sort data=sim_list;
  by transcript_id;
run;

data event_summary;
  merge ea_list (in=in1) sim_list (in=in2);
  by transcript_id;
  if in1 then flag_in_reduced_ref=1; else flag_in_reduced_ref=0;
  if in2 then flag_in_simulated_list=1; else flag_in_simulated_list=0;
run;

proc freq data=event_summary noprint;
  tables flag_in_simulated_list*flag_in_reduced_ref / out=xs_count;
run;

proc print data=xs_count;
run;

%mend;

%countEA(1,0);
%countEA(0.75,0);
%countEA(0.75,5);


/*
 flag_in_     flag_in_
simulated_    reduced_
   list          ref      COUNT

     0            1        1048
     1            0         325
     1            1         142


142 simulated transcripts kept
325 simulated transcripts missed
1048 non-simulated transcripts added

75% APN0
  flag_in_     flag_in_
 simulated_    reduced_
    list          ref      COUNT

      0            1        3388
      1            0         165
      1            1         302
302 simulated transcripts kept
165 simulated transcripts missed
3388 non-simulated transcripts added


75% APN5
  flag_in_     flag_in_
 simulated_    reduced_
    list          ref      COUNT

      0            1        2557
      1            0         170
      1            1         297
297 simulated transcripts kept
170 simulated transcripts missed
2557 non-simulated transcripts added

*/
