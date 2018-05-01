ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";


/* For each set of transcripts calculate the proportion of transcripts with a diff_log_tpm of 2 */


%macro propdiff(datain);

data tpm_data;
   set event.rsem_&datain.;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id  mean_tpm tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
run;

data calc_diff;
  retain transcript_id tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 tpm_bin mean_tpm;
  set tpm_data;
  diff_tpm=(tpm_nsc1-tpm_nsc2);
  mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
  diff_log_tpm=(log_tpm_nsc1-log_tpm_nsc2);
  if abs(diff_log_tpm) ge 1.5 then flag_diff_ge15=1;
  else flag_diff_ge15=0;
run;

proc freq data=calc_diff;
   tables flag_diff_ge15;
run;


*proc freq data=calc_diff;
*   where mean_tpm>0;
*   tables flag_diff_ge15;
*run;
%mend;

%propdiff(pacbio_all);
%propdiff(refseq_all);
%propdiff(events_exp_any);
%propdiff(events_exp_100perc);
%propdiff(events_exp_75perc);
%propdiff(events_exp_75perc_apn5);
%propdiff(events_exp_75perc_apn10);
%propdiff(events_exp_75perc_apn10);
%propdiff(events_exp_100perc_apn5);
%propdiff(events_exp_100perc_apn10);
%propdiff(events_exp_50perc);
%propdiff(events_exp_50perc_apn5);
%propdiff(events_exp_50perc_apn10);
%propdiff(pacbio_apn0);
%propdiff(pacbio_apn5);
%propdiff(pacbio_apn10);
/*
		ALL			
		N(all)	N(diff)	%(diff)	
PacBio All	16104	254	1.58
RefSeq All	128631	329	0.26
Events Any	73535	266	0.36
EA 100%		14734	149	1.01
EA 75%		34622	259	0.75
EA 75%, APN5	13557	183	1.33
EA 75%, APN10	9879	149	1.51
EA 100%, APN5	3815	56	1.47
EA 100%, APN10	2606	41	1.57
EA 50%		45883	269	0.59
EA 50%, APN5	20336	227	1.12
EA 50%, APN10	15102	181	1.20
PacBio APN0	6286	63	1.00
PacBio APN5	6078	63	1.04
PacBio APN10	5884	64	1.09
*/



