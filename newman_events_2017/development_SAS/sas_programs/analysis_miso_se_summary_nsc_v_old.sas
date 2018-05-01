ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Summarize MISO skipped exons */

data miso_se;
   set event.miso_bin_diffs_and_bayes;
   if NSC_OLD_diff ne .;
   keep event_name NSC_OLD_diff NSC_OLD_Bayes_factor;
run;

data miso_diff_se;
   set miso_se;
   if NSC_OLD_diff ge 0.2 and NSC_OLD_Bayes_factor ge 10 then flag_miso_diff_se_bf10=1;
   else flag_miso_diff_se_bf10=0;
   if NSC_OLD_diff ge 0.2 and NSC_OLD_Bayes_factor ge 5  then flag_miso_diff_se_bf5=1;
   else flag_miso_diff_se_bf5=0;
run;

data miso2gene;
  set event.miso_se2gene;
run;

proc sort data=miso_diff_se;
   by event_name;
proc sort data=miso2gene;
   by event_name;
run;

data miso2gene_w_diff;
  merge miso2gene (in=in1) miso_diff_se (in=in2);
  by event_name;
  if in2;
run;

proc sort data=miso2gene_w_diff;
  by ens_gene_id;
proc means data=miso2gene_w_diff noprint;
  by ens_gene_id;
  var flag_miso_diff_se_bf10 flag_miso_diff_se_bf5;
  output out=miso_gene_w_diff_se_cnt sum(flag_miso_diff_se_bf10)=num_diff_se_bf10 sum(flag_miso_diff_se_bf5)=num_diff_se_bf5;
run;


data ens2refseq;
   set event.ensembl2refseq_gene_id;
run;

proc sort data=ens2refseq;
  by ens_gene_id;
proc sort data=miso_gene_w_diff_se_cnt;
  by ens_gene_id;
run;

data miso_diff_se_cnt2refseq;
  merge miso_gene_w_diff_se_cnt (in=in1) ens2refseq (in=in2);
  by ens_gene_id;
  if in1 and in2 then do;
         flag_refseq_match=1;
         output; end;
  else if in1 then do;
         flag_refseq_match=0;
         output; end;
run;

data event.miso_diff_se_refseq;
  set miso_diff_se_cnt2refseq;
  drop _TYPE_;
  rename _FREQ_ = num_se_events_total;
run;

