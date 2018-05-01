ods listing; ods html close;
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/*
Estimate APN per transcript from event APN and compare to TPM estimates
do this for the case-only (all transcripts and all reduced) then subset the all reduced
for the genes with kappa <0.6

Do this for fragments, since I don't think I still have the splicing counts

(1) Frag APN * length
(2) Sum frags per transcript for length and depth
(3) Recalc APN


For transcripts in the reduced transcriptome, I want to estimate the APN from fragment coverage and then compare the TPM estimate with the APN.

I will do this for all transcripts in the reduced set, then only the MEI(s) */

/* Estimate APN */

data frag_counts;
  set eventloc.t1d_case_only_fragment_counts;
  keep name fragment_id cell_type apn;
  rename fragment_id=fragment_coord;
run;

data frag_id;
  length fragment_coord $30.;
  set hg19.hg19_aceview_exon_fragment_info;
  fragment_coord=catx(":",chr,fragment_start,fragment_end);
  fragment_length=fragment_start-fragment_end;
  keep fragment_id fragment_coord transcript_id fragment_length;
run;

proc sort data=frag_counts;
  by fragment_coord;
proc sort data=frag_id;
  by fragment_coord;
run;

data frag_counts_w_xs;
  merge frag_id (in=in1) frag_counts (in=in2);
  by fragment_coord;
  if in1 and in2;
run;

/* Expand transcripts */

data stack_xs;
  length transcript_id2 $60.;
  set frag_counts_w_xs;
  fragment_length=abs(fragment_length);
  do i=1 by 1 while(scan(transcript_id,i,"|") ^= "" );
      transcript_id2=scan(transcript_id,i,"|");
      output; end;
  keep name cell_type apn fragment_id fragment_length transcript_id2;
  rename transcript_id2=transcript_id;
run;

/* Get list of transcripts to estimate */

data xs2keep;
   set eventloc.hg19_rsem_75perc_apn5_xscripts;
   keep transcript_id;
run;

proc sort data=xs2keep nodup;
  by transcript_id;
proc sort data=stack_xs;
  by transcript_id;
run;

data stack_xs_subset;
  merge xs2keep (in=in1) stack_xs (in=in2);
  by transcript_id;
  if in1 and in2;
run;

data calc_depth;
  set stack_xs_subset;
  depth=apn*fragment_length;
run;

proc sort data=calc_depth;
  by name cell_type transcript_id fragment_id;
proc means data=calc_depth noprint;
  by name cell_type transcript_id;
  var fragment_length depth;
  output out=sum_frags_by_xs sum=;
run;

data recalc_apn;
  set sum_frags_by_xs;
  apn=depth/fragment_length;
  log_apn=log(apn+1);
  drop _TYPE_ _FREQ_;
run;

/* Make permenant */

data eventloc.hg19_rsem_reduced_xs_apn_est;
   set recalc_apn;
run;


