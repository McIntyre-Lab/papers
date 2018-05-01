ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';

/* Import eXpress results. The results format has the following columns:
    bundle_id          ID of bundle the target belongs to. A bundle is defined
                       as the transitive closure of targets that share multi-mapping reads.
    target_id          transcript ID
    length             length of transcript sequence
    eff_length         length of transcript seq, adjusted for fragment bias
    tot_counts         total number of fragments mapping
    uniq_counts        number of fragments uniquely mapping
    est_counts         estimated number of fragments generated from transcript
    eff_counts         est_counts, adjusted for fragment and length bias
    ambig_distr_alpha  Alpha parameter for posterior beta-binomial distribution fit to ambig reads
    ambig_distr_beta   Beta parameter for posterior beta-binomial distribution fit to ambig reads
    fpkm               Fragments per kilobase per million mapped
    fpkm_conf_low      95% lower FPKM
    fpkm_conf_high     95% upper FPKM
    solvable           Binary indicator for whether likelihood function has a unique maximum
    tpm                transcripts per million
 */

%macro processXprs(xsList,outName);

   %macro import_data(sample);
   proc import datafile="!MCLAB/event_analysis/analysis_output/express_output/&sample._&xsList./results.xprs"
        out=event.xprs_&sample._&outName. dbms=tab replace;
        guessingrows=130000;
   run;
  
   data xprs_&sample._2;
      length sample_id $4.;
      set event.xprs_&sample._&outName. ;
      sample_id="&sample.";
      log_tpm=log(TPM + 1);
      /* Bin transcripts */
      if log_tpm=0 then tpm_bin=0;
      else if log_tpm < 0.5 then tpm_bin=1;
      else if log_tpm < 2 then tpm_bin=2;
      else if log_tpm < 4 then tpm_bin=3;
      else tpm_bin=4;
      keep sample_id target_id log_tpm tpm tpm_bin;
      rename target_id=transcript_id;
   run;     

   %mend;
   %import_data(NSC1);
   %import_data(NSC2);

   data stack_xprs;
      set xprs_NSC1_2 xprs_NSC2_2;
   run;

   /* Calc CV by sample and bin */
   proc sort data=stack_xprs;
      by sample_id tpm_bin;
   proc means data=stack_xprs noprint;
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
      rename NSC1=NSC1_CV NSC2=NSC2_CV;
      drop _NAME_ ;
   run;

   data count_sbys2;
      set count_sbys;
      rename NSC1=NSC1_count NSC2=NSC2_count;
      drop _NAME_ ;
   run;

   data stats_express_&outName.;
      length xscript_set $32.;
      merge count_sbys2 cv_sbys2;
      by tpm_bin;
      xscript_set="&xsList.";
   run;

   /* Make permenant */
   data event.stats_express_&outName.;
     set stats_express_&outName.;
   run;
%mend;

%processXprs(refseq_mm10,refseq_all);
%processXprs(refseq_mm10_exp_transcripts_100perc_dtct,events_100prc_apn0);
%processXprs(refseq_mm10_exp_transcripts_75perc_dtct_apn5,events_75prc_apn5);


