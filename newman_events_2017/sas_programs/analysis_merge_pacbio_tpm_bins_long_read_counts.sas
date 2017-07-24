ods listing; ods html close;
libname event "!MCLAB/event_analysis/sas_data";

/* Merge PacBio long read counts with TPM bins for each set of PB transcripts used */

data lr_counts;
   set event.pacbio_lr_read_counts;
   mean_npc_fl=mean(reads_FL_NPC1,reads_FL_NPC2);
   mean_npc_all=mean((reads_FL_NPC1+reads_nonFL_NPC1),(reads_FL_NPC2+reads_nonFL_NPC2));
   keep pacbio_id mean_npc_fl mean_npc_all;
   rename pacbio_id=transcript_id;
run;


%macro bins(datain);

data tpm_data;
   set event.rsem_&datain.;
   mean_tpm=(tpm_nsc1+tpm_nsc2)/2;
   keep transcript_id mean_tpm;
run;

proc means data=tpm_data noprint;
  where mean_tpm > 0;
  var mean_tpm;
  output out=tpm_distrib min=min q1=q1 q3=q3 max=max;
run;

%local tpmMin tpmQ1 tpmQ3 tpmMax;

data _null_;
   set tpm_distrib;
   call symputx('tpmMin',min);
   call symputx('tpmQ1',q1);
   call symputx('tpmQ3',q3);
   call symputx('tpmMax',max);
run;

data bins_&datain.;
  retain transcript_id tpm_bin_&datain.;
  set tpm_data;
  if mean_tpm le &tpmQ1. then tpm_bin_&datain.=1;
  else if mean_tpm le &tpmQ3. then tpm_bin_&datain.=2;
  else tpm_bin_&datain.=3;
  keep transcript_id tpm_bin_&datain.;
run;

proc sort data=bins_&datain.;
   by transcript_id;
run;

%mend;

%bins(pacbio_all);
%bins(pacbio_apn0);
%bins(pacbio_apn5);
%bins(pacbio_apn10);

proc sort data=lr_counts;
   by transcript_id;
run;

data lr_w_bins;
  merge lr_counts (in=in1) bins_pacbio_all (in=in2) bins_pacbio_apn0 (in=in3) 
        bins_pacbio_apn5 (in=in4) bins_pacbio_apn10 (in=in5);
  by transcript_id;
  if not in2 then tpm_bin_pacbio_all=0;
  if not in3 then tpm_bin_pacbio_apn0=0;
  if not in4 then tpm_bin_pacbio_apn5=0;
  if not in5 then tpm_bin_pacbio_apn10=0;
run;


/* Export for plots */

proc export data=lr_w_bins
     outfile="!MCLAB/event_analysis/analysis_output/pacbio_lr_count_w_tpm_bin.csv"
     dbms=csv replace;
run;

