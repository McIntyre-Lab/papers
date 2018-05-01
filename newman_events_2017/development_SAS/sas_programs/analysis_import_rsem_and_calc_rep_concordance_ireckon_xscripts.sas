ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";

/* Macro to import and save RSEM results, and calculate correlation between replicates */



/* Import data */
proc import datafile="!MCLAB/event_analysis/alignment_output/rsem_ireckon_transcripts_npc/NSC1.isoforms.results" out=rsem_NSC1 dbms=tab replace;
guessingrows=24000;
run;

proc import datafile="!MCLAB/event_analysis/alignment_output/rsem_ireckon_transcripts_npc/NSC2.isoforms.results" out=rsem_NSC2 dbms=tab replace;
guessingrows=24000;
run;



data rsem_NSC1_2;
  set rsem_NSC1;
  keep transcript_id tpm posterior_mean_count;
  rename tpm=tpm_nsc1
         posterior_mean_count=post_mean_nsc1;
run;

data rsem_NSC2_2;
  set rsem_NSC2;
  keep transcript_id tpm posterior_mean_count;
  rename tpm=tpm_nsc2
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
            post_mean_nsc1=0;
            end;
  if in2 then flag_in_rep2=1;
  	     else do;
            flag_in_rep2=0;
            tpm_nsc2=0;
            post_mean_nsc2=0;
            end;
run;

data rsem_NSC_2;
   set rsem_NSC;
   log_tpm_nsc1=log(tpm_nsc1+1);
   log_tpm_nsc2=log(tpm_nsc2+1);

run;

proc corr data=rsem_NSC_2 pearson;
  var log_tpm_nsc1 log_tpm_nsc2;
run;

/*

  Pearson Correlation Coefficients, N = 23738
           Prob > |r| under H0: Rho=0

                      log_tpm_      log_tpm_
                          nsc1          nsc2

    log_tpm_nsc1       1.00000       0.96387
                                      <.0001

    log_tpm_nsc2       0.96387       1.00000
                        <.0001

*/

/* Make permenant */

data event.rsem_ireckon_transcripts_npc;
    set rsem_NSC_2;
run;

/* For transcripts by type (novel, IR, known Refseq, unspliced), what is the TPM range? */

data rsem_NSC;
  set event.rsem_ireckon_transcripts_npc;
  length transcript_type $100.;
  length sample_source $4.;
  length gene_id $100.;
  length transcript_id2 $100.;
  sample_source=scan(transcript_id,1,"_");
  gene_id=scan(transcript_id,2,"_");
  if index(transcript_id,"NM_") > 0 or index(transcript_id,"NR_") > 0 or
     index(transcript_id,"XM_") > 0 or index(transcript_id,"XR_") > 0 
     then transcript_id2=catx("_",scan(transcript_id,3,"_"),scan(transcript_id,4,"_"));
  else transcript_id2=catx("_",scan(transcript_id,3,"_"));

  if index(transcript_id2,"unspliced") > 0 then transcript_type="unspliced";
  else if index(transcript_id2,"novel") > 0 then transcript_type="novel";
  else if index(transcript_id2,"Intron") > 0 then transcript_type="intron_retention";
  else if index(transcript_id,"NM_") > 0 or index(transcript_id,"NR_") > 0 or
     index(transcript_id,"XM_") > 0 or index(transcript_id,"XR_") > 0 
     then transcript_type="Known refseq";
  else transcript_type="unknown";
  mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
  mean_log_tpm=(log_tpm_nsc1+log_tpm_nsc2)/2;
run;

proc freq data=rsem_NSC;
  table transcript_type;
run;

/*
   transcript_type     Frequency
   ------------------------------
   Known refseq            1519
   intron_retention         109
   novel                  16879
   unspliced               5231


*/

proc sort data=rsem_NSC;
  by transcript_type;
proc means data=rsem_nsc noprint;
  by transcript_type;
  var tpm_nsc1 tpm_nsc2 log_tpm_nsc1 log_tpm_nsc2 mean_tpm mean_log_tpm;
  output out=rsem_ir_stats;
run;

proc print data=rsem_ir_stats;
run;


/* For known transcripts, how well do the TPM estimates compare between EA and iReckon ? 

Do for the 100% APN0 and 75% APN5 sets
*/

data rsem_xs100;
  set event.rsem_events_exp_100perc;
  mean_log_tpm_ea=(log(tpm_nsc1+1) + log(tpm_nsc2+1)) / 2;
  log_tpm_nsc1_ea=log(tpm_nsc1+1);
  log_tpm_nsc2_ea=log(tpm_nsc2+1);
  keep transcript_id mean_log_tpm_ea log_tpm_nsc1_ea log_tpm_nsc2_ea;
  rename transcript_id=transcript_id2;
run;

data rsem_xs75;
  set event.rsem_events_exp_75perc_apn5;
  mean_log_tpm_ea=(log(tpm_nsc1+1) + log(tpm_nsc2+1)) / 2;
  log_tpm_nsc1_ea=log(tpm_nsc1+1);
  log_tpm_nsc2_ea=log(tpm_nsc2+1);
  keep transcript_id mean_log_tpm_ea log_tpm_nsc1_ea log_tpm_nsc2_ea;
  rename transcript_id=transcript_id2;
run;

proc sort data=rsem_xs100;
  by transcript_id2;
proc sort data=rsem_xs75;
  by transcript_id2;
proc sort data=rsem_nsc;
  by transcript_id2;
run;

data rsem_ir_v_xs100;
  merge rsem_nsc (in=in1) rsem_xs100 (in=in2);
  by transcript_id2;
  if in1 and in2; *145 transcripts in common;
run;



data rsem_ir_v_xs75;
  merge rsem_nsc (in=in1) rsem_xs75 (in=in2);
  by transcript_id2;
  if in1 and in2; *161 transcripts in common;
run;

ods graphics on;
proc sgplot data=rsem_ir_v_xs100;
   scatter x=mean_log_tpm_ea y=mean_log_tpm;
run;

proc sgplot data=rsem_ir_v_xs75;
   scatter x=mean_log_tpm_ea y=mean_log_tpm;
run;

proc corr data=rsem_ir_v_xs100 pearson;
  var mean_log_tpm_ea mean_log_tpm;
run;

proc corr data=rsem_ir_v_xs75 pearson;
  var mean_log_tpm_ea mean_log_tpm;
run;

/*
                          mean_
                       log_tpm_         mean_
                             ea       log_tpm

  mean_log_tpm_ea       1.00000       0.58453
                                       <.0001

  mean_log_tpm          0.58453       1.00000
                         <.0001



                          mean_
                       log_tpm_         mean_
                             ea       log_tpm

  mean_log_tpm_ea       1.00000       0.57437
                                       <.0001

  mean_log_tpm          0.57437       1.00000
                         <.0001




*/

