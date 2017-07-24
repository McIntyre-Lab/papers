
/* Import MISO comparisons. Put all these together so LMM and I can look at them and decide what to do */

ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';


/* Import data */

%macro import_miso(sample1,sample2,counter);

    data WORK.MISO_&sample1._&sample2.    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile
"!MCLAB/event_analysis/analysis_output/MISO/comparisons/miso_output_SE_&sample1._vs_miso_output_SE_
&sample2./bayes-factors/miso_output_SE_&sample1._vs_miso_output_SE_&sample2..miso_bf" delimiter='09'x
MISSOVER
 DSD lrecl=32767 firstobs=2 ;
        informat event_name $83. ;
        informat sample1_posterior_mean best32. ;
        informat sample1_ci_low best32. ;
        informat sample1_ci_high best32. ;
        informat sample2_posterior_mean best32. ;
        informat sample2_ci_low best32. ;
        informat sample2_ci_high best32. ;
        informat diff best32. ;
        informat bayes_factor best32. ;
        informat isoforms $448. ;
        informat sample1_counts $34. ;
        informat sample1_assigned_counts $15. ;
        informat sample2_counts $34. ;
        informat sample2_assigned_counts $15. ;
        informat chrom $5. ;
        informat strand $1. ;
        informat mRNA_starts $20. ;
        informat mRNA_ends $20. ;
        format event_name $83. ;
        format sample1_posterior_mean best12. ;
        format sample1_ci_low best12. ;
        format sample1_ci_high best12. ;
        format sample2_posterior_mean best12. ;
        format sample2_ci_low best12. ;
        format sample2_ci_high best12. ;
        format diff best12. ;
        format bayes_factor best12. ;
        format isoforms $448. ;
        format sample1_counts $34. ;
        format sample1_assigned_counts $15. ;
        format sample2_counts $34. ;
        format sample2_assigned_counts $15. ;
        format chrom $5. ;
        format strand $1. ;
        format mRNA_starts $20. ;
        format mRNA_ends $20. ;
   input
               event_name $
               sample1_posterior_mean
               sample1_ci_low
               sample1_ci_high
               sample2_posterior_mean
               sample2_ci_low
               sample2_ci_high
               diff
               bayes_factor
               isoforms $
               sample1_counts $
               sample1_assigned_counts $
               sample2_counts $
               sample2_assigned_counts $
               chrom $
               strand $
               mRNA_starts $
               mRNA_ends $
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

data MISO_&sample1._&sample2._2;
  set MISO_&sample1._&sample2.;
  rename sample1_posterior_mean=&sample1._posterior_mean_&counter.
         sample1_ci_low=&sample1._ci_low_&counter.
         sample1_ci_high=&sample1._ci_high_&counter.
         sample2_posterior_mean=&sample2._posterior_mean_&counter.
         sample2_ci_low=&sample2._ci_low_&counter.
         sample2_ci_high=&sample2._ci_high_&counter.
         diff=&sample1._&sample2._diff_&counter.
         bayes_factor=&sample1._&sample2._bayes_factor_&counter.
         isoforms=&sample1._&sample2._isoforms_&counter.
         sample1_counts=&sample1._counts_&counter.
         sample1_assigned_counts=&sample1._assigned_counts_&counter.
         sample2_counts=&sample2._counts_&counter.
         sample2_assigned_counts=&sample2._assigned_counts_&counter.
         chrom=&sample1._&sample2._chrom_&counter.
         strand=&sample1._&sample2._strand_&counter.
         mRNA_starts=&sample1._&sample2._mRNA_starts_&counter.
         mRNA_ends=&sample1._&sample2._mRNA_ends_&counter.;
run;

%mend;


%import_miso(NSC1,NSC2,1);
%import_miso(NSC1,OLD1,2);
%import_miso(NSC1,OLD2,3);
%import_miso(NSC2,OLD1,4);
%import_miso(NSC2,OLD2,5);
%import_miso(OLD1,OLD2,6);
%import_miso(NSC,OLD,7);

proc sort data=MISO_NSC1_NSC2_2;
   by event_name;
proc sort data=MISO_NSC1_OLD1_2;
   by event_name;
proc sort data=MISO_NSC1_OLD2_2;
   by event_name;
proc sort data=MISO_NSC2_OLD1_2;
   by event_name;
proc sort data=MISO_NSC2_OLD2_2;
   by event_name;
proc sort data=MISO_OLD1_OLD2_2;
   by event_name;
proc sort data=MISO_NSC_OLD_2;
   by event_name;
run;

data miso_all_results;
  merge MISO_NSC1_NSC2_2 MISO_NSC1_OLD1_2 MISO_NSC1_OLD2_2
        MISO_NSC2_OLD1_2 MISO_NSC2_OLD2_2 MISO_OLD1_OLD2_2
        MISO_NSC_OLD_2 ;
  by event_name;
run;


/* Make permenant */

data event.miso_all_results_nsc_v_old;
   set miso_all_results;
run;


