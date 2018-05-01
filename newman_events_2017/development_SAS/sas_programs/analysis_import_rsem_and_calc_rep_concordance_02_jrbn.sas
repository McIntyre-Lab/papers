ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";

/* Macro to import and save RSEM results, and calculate correlation between replicates */

%macro rsem(datain,outname);


/* Import data */
proc import datafile="!MCLAB/event_analysis/analysis_output/rsem_output/&datain._NSC1.isoforms.results" out=rsem_NSC1 dbms=tab replace;
guessingrows=128632;
run;

proc import datafile="!MCLAB/event_analysis/analysis_output/rsem_output/&datain._NSC2.isoforms.results" out=rsem_NSC2 dbms=tab replace;
guessingrows=128632;
run;



data rsem_NSC1_2;
  set rsem_NSC1;
  keep transcript_id tpm fpkm posterior_mean_count;
  rename tpm=tpm_nsc1
         fpkm=fpkm_nsc1
         posterior_mean_count=post_mean_nsc1;
run;

data rsem_NSC2_2;
  set rsem_NSC2;
  keep transcript_id tpm fpkm posterior_mean_count;
  rename tpm=tpm_nsc2
         fpkm=fpkm_nsc2
         posterior_mean_count=post_mean_nsc2;
run;

proc sort data=rsem_NSC1_2;
   by transcript_id;
proc sort data=rsem_NSC2_2;
   by transcript_id;
run;

data rsem_NSC;
  merge rsem_NSC1_2 (in=in1) rsem_NSC2_2 (in=in2);
  by transcript_id;
  if in1 then flag_in_rep1=1;
         else do;
            flag_in_rep1=0;
            tpm_nsc1=0;
            fpkm_nsc1=0;
            post_mean_nsc1=0;
            end;
  if in2 then flag_in_rep2=1;
  	     else do;
            flag_in_rep2=0;
            tpm_nsc2=0;
            fpkm_nsc2=0;
            post_mean_nsc2=0;
            end;
run;

data rsem_NSC_2;
   set rsem_NSC;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   log_fpkm_nsc1=log(fpkm_nsc1+1);
   log_fpkm_nsc2=log(fpkm_nsc2+1);
run;

proc corr data=rsem_NSC pearson;
  var fpkm_nsc1 fpkm_nsc2 log_fpkm_nsc1 log_fpkm_nsc2;
run;

proc corr data=rsem_NSC spearman;
  var fpkm_nsc1 fpkm_nsc2 log_fpkm_nsc1 log_fpkm_nsc2;
run;

/* Make permenant */

data event.rsem_&outname.;
    set rsem_NSC;
run;

data rsem_NSC;
set event.rsem_&outname.;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);
   log_fpkm_nsc1=log(fpkm_nsc1+1);
   log_fpkm_nsc2=log(fpkm_nsc2+1);
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
run;


proc corr data=rsem_NSC pearson;
  var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
run;

proc corr data=rsem_NSC spearman;
  var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2;
run;

proc print data=rsem_NSC(keep=transcript_id mean_tpm);
run;



%mend;


%rsem(refseq_mm10,refseq_all);
%rsem(pacbio_mm10,pacbio_all);
%rsem(refseq_mm10_exp_transcripts_100perc_dtct,events_exp_100perc);
%rsem(refseq_mm10_exp_transcripts_75perc_dtct_apn5,events_exp_75perc_apn5);
%rsem(refseq_mm10_exp_transcripts,events_exp_any);
%rsem(refseq_mm10_exp_transcripts_75perc_dtct,events_exp_75perc);
%rsem(refseq_mm10_exp_transcripts_75perc_dtct_apn10,events_exp_75perc_apn10);
%rsem(refseq_mm10_exp_transcripts_100perc_dtct_apn5,events_exp_100perc_apn5);
%rsem(refseq_mm10_exp_transcripts_100perc_dtct_apn10,events_exp_100perc_apn10);
%rsem(refseq_mm10_exp_transcripts_50perc_dtct,events_exp_50perc);
%rsem(refseq_mm10_exp_transcripts_50perc_dtct_apn5,events_exp_50perc_apn5);
%rsem(refseq_mm10_exp_transcripts_50perc_dtct_apn10,events_exp_50perc_apn10);
%rsem(pacbio_mm10_apn0,pacbio_apn0);
%rsem(pacbio_mm10_apn5,pacbio_apn5);
%rsem(pacbio_mm10_apn10,pacbio_apn10);

%rsem(cav1_set1,cav1_set1);
%rsem(cav1_set2,cav1_set2);
%rsem(cav1_set3,cav1_set3);

