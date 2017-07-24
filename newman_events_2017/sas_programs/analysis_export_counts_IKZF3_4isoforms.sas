ods listing; ods html close;

/* Prepare data for two figures:

(1) Mean TPM by cell type and IKZF3 transcript

(2) TPM by cell type, transcript and subject
 should look like:
Fig 1:
TranscriptID, CD4, CD8, CD19
IKZF3.#Aug10, mean values

Fig 2:

SubjectID	CD4	CD8	CD19
		xs1,xs2,xs3

*/

data ikzf3_tpm;
  set eventloc.hg19_rsem_75perc_apn5_xscripts;
  where transcript_id="IKZF3.aAug10"
     or transcript_id="IKZF3.bAug10"
     or transcript_id="IKZF3.iAug10"
     or transcript_id="IKZF3.oAug10";
  log_tpm=log(tpm+1);
run;

data design;
   set con.design_by_subject_new;
   /* Drop these samples -- low coverage */
   if name =  '2009-PC-0221' then delete; *sample 75 cd8;
   if name =  '2009-PC-0144' then delete; *sample 48 cd4;
   if name =  '2009-PC-0236' then delete; *sample 80;
   if name =  '2009-PC-0237' then delete; *sample 80;
   if name =  '2009-PC-0235' then delete; *sample 80; 
   keep subject_id cell_type library;
run;

proc sort data=ikzf3_tpm;
   by library;
proc sort data=design;
   by library;
run;

data ikzf3_w_key;
  merge design (in=in1) ikzf3_tpm (in=in2);
  by library;
  if in1 and in2;
run;

/* Calculate mean by cell type, transcript */

proc sort data=ikzf3_w_key;
   by transcript_id cell_type;
proc means data=ikzf3_w_key noprint;
  by transcript_id cell_type;
  var tpm;
  output out=mean_tpm_by_xs_cell mean=;
run;

data mean_tpm_by_xs_cell2;
  set mean_tpm_by_xs_cell;
  log_tpm=log(tpm+1);
run;

proc sort data=mean_tpm_by_xs_cell2;
  by transcript_id cell_type;
proc transpose data=mean_tpm_by_xs_cell2 out=mean_tpm_sbys;
  by transcript_id;
  id cell_type;
  var log_tpm;
run;
 
proc export data=mean_tpm_sbys outfile="!MCLAB/event_analysis/analysis_output/IKZF3_mean_tpm_by_cell_transcript.csv"
   dbms=csv replace;
run;

/* Now put by-subject estimates SbyS */

data ikzf3_new_header;
  length cell_by_xs $20.;
  set ikzf3_w_key;
  cell_by_xs=catx("_",cell_type,scan(transcript_id,2,"."));
run;

proc sort data=ikzf3_new_header;
   by subject_id cell_by_xs;
proc transpose data=ikzf3_new_header out=ikzf3_by_subj_sbys;
   by subject_id;
   id cell_by_xs;
   var log_tpm;
run;

/* export */

proc export data=ikzf3_by_subj_sbys outfile="!MCLAB/event_analysis/analysis_output/IKZF3_tpm_by_subject_cell_transcript.csv"
    dbms=csv replace;
run;

