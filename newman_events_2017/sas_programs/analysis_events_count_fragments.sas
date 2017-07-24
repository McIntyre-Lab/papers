ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count the number of exon fragments that are expressed at APN>0 */

data fragments;
  set event.flag_fragment_on;
  keep fragment_id flag_fragment_nsc_on;
  rename flag_fragment_nsc_on=flag_fragment_on;
run;


proc freq data=fragments;
   tables flag_fragment_on;
run;


/*

                                                 Cumulative    Cumulative
    flag_fragment_on    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0      141755       39.35        141755        39.35
                   1      218481       60.65        360236       100.00

*/


data event.fragments_on_apn_gt0;
   set fragments;
run;

