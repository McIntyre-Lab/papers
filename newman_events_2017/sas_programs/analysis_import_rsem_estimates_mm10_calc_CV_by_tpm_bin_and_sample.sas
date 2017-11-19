ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";

/* Macro to import and save RSEM results, and calculate CV per transcript bin */

%macro rsem(datain,outname);
/* Import data */
proc import datafile="!MCLAB/event_analysis/analysis_output/old/rsem_output/&datain._NSC1.isoforms.results" out=rsem_NSC1 dbms=tab replace;
guessingrows=128632;
run;
proc import datafile="!MCLAB/event_analysis/analysis_output/old/rsem_output/&datain._NSC2.isoforms.results" out=rsem_NSC2 dbms=tab replace;
guessingrows=128632;
run;
proc import datafile="!MCLAB/event_analysis/analysis_output/old/rsem_output/&datain._OLD1.isoforms.results" out=rsem_OLD1 dbms=tab replace;
guessingrows=128632;
run;
proc import datafile="!MCLAB/event_analysis/analysis_output/old/rsem_output/&datain._OLD2.isoforms.results" out=rsem_OLD2 dbms=tab replace;
guessingrows=128632;
run;

/* Stack data */

%macro bin_xscripts(sample);
data rsem_&sample._2;
  length sample_id $4.;
  set rsem_&sample.;
  sample_id="&sample.";
  log_tpm=log(TPM + 1);
  /* Bin transcripts */
  if log_tpm=0 then tpm_bin=0;
  else if log_tpm < 0.5 then tpm_bin=1;
  else if log_tpm < 2 then tpm_bin=2;
  else if log_tpm < 4 then tpm_bin=3;
  else tpm_bin=4;
  keep sample_id transcript_id log_tpm tpm tpm_bin;
run;
%mend;

%bin_xscripts(NSC1);
%bin_xscripts(NSC2);
%bin_xscripts(OLD1);
%bin_xscripts(OLD2);

data stack_rsem;
   set rsem_NSC1_2 rsem_NSC2_2 rsem_OLD1_2 rsem_OLD2_2;
run;

/* Calc CV by sample and bin */
proc sort data=stack_rsem;
  by sample_id tpm_bin;
proc means data=stack_rsem noprint;
  by sample_id tpm_bin;
  var log_tpm;
  output out=varstats_by_bin_sample cv=tpm_cv;
run;

/* Now I need to transpose this, as I want bin as rows and samples as columns */

proc sort data=varstats_by_bin_sample;
   by tpm_bin sample_id;
proc transpose data=varstats_by_bin_sample out=cv_sbys;
  by tpm_bin;
  id sample_id;
  var tpm_cv;
run;

proc transpose data=varstats_by_bin_sample out=count_sbys;
  by tpm_bin;
  id sample_id;
  var _FREQ_;
run;

data cv_sbys2;
   set cv_sbys;
   rename NSC1=NSC1_CV NSC2=NSC2_CV OLD1=OLD1_CV OLD2=OLD2_CV;
   drop _NAME_ ;
run;

data count_sbys2;
   set count_sbys;
   rename NSC1=NSC1_count NSC2=NSC2_count OLD1=OLD1_count OLD2=OLD2_count;
   drop _NAME_ ;
run;

data stats_&outname.;
   length xscript_set $32.;
   merge count_sbys2 cv_sbys2;
   by tpm_bin;
   xscript_set="&outname.";
run;

%mend;


%rsem(pacbio_mm10,pacbio_all);
%rsem(refseq_mm10,refseq_all);
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


data all_stats_stacked;
   set stats_: ;
run;

proc export data=all_stats_stacked
     outfile="!MCLAB/event_analysis/analysis_output/CV_and_counts_by_tpm_bin_and_sample_alldata.csv"
     dbms=csv replace;
run;


   
