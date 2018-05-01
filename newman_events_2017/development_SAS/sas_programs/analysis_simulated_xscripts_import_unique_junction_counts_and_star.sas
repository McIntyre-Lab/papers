/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the benchmarked simulated data, import the STAR junction results and EA junction coverage counts and check
   to see which method captures the junctions of these genes */

/* Import STAR */

%macro importSTAR(sample);

proc import datafile="!MCLAB/event_analysis/alignment_output/aln_mm10_star_simulation_100bp/&sample.SJ.out.tab"
    out=&sample._star dbms=tab replace;
    getnames=no; guessingrows=max;
run;

data &sample._star2;
  length sample_id $15.;
  set &sample._star;
  sample_id="&sample.";
  rename VAR1=chr
         VAR2=intron_start
         VAR3=intron_stop
         VAR4=strand
         VAR5=intron_motif_type
         VAR6=flag_junction_annotated
         VAR7=num_unique_mapped_reads
         VAR8=num_multimapped_reads
         VAR9=max_overhang
         ;
run;

%mend;

%importSTAR(sim1_test1);
%importSTAR(sim1_test2);
%importSTAR(sim2_test1);
%importSTAR(sim2_test2);
%importSTAR(sim3_test1);
%importSTAR(sim3_test2);

/* Import coverage counts */

%macro importCC(sample);

    data &sample._counts    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "!MCLAB/event_analysis/alignment_output/coverage_counts_splicing_simul/cvrg_cnts_&sample..csv" delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $10. ;
       informat event_id $16. ;
       informat mapped_reads best32. ;
       informat read_length best32. ;
       informat region_length best32. ;
       informat region_depth best32. ;
       informat reads_in_region best32. ;
       informat apn best32. ;
       informat rpkm best32. ;
       informat mean best32. ;
       informat std best32. ;
       informat cv best32. ;
       format sample_id $10. ;
       format event_id $16. ;
       format mapped_reads best12. ;
       format read_length best12. ;
       format region_length best12. ;
       format region_depth best12. ;
       format reads_in_region best12. ;
       format apn best12. ;
       format rpkm best12. ;
       format mean best12. ;
       format std best12. ;
       format cv best12. ;
    input
               sample_id $
               event_id $
               mapped_reads
               read_length
               region_length
               region_depth
               reads_in_region
               apn
               rpkm
               mean
               std
               cv
   ;
   if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
   run;

%mend;

%importCC(sim1_test1);
%importCC(sim1_test2);
%importCC(sim2_test1);
%importCC(sim2_test2);
%importCC(sim3_test1);
%importCC(sim3_test2);

/* Stack outputs and save */

data star_juncs;
   set sim1_test1_star2 sim1_test2_star2 sim2_test1_star2
       sim2_test2_star2 sim3_test1_star2 sim3_test2_star2;
run;

data junc_counts;
   set sim1_test1_counts sim1_test2_counts
       sim2_test1_counts sim2_test2_counts
       sim3_test1_counts sim3_test2_counts;
run;

/* For catalog junctions, I really only want to keep the junctions from the transcripts selected, but I also want
   to double check that reads didn't go to the wrong junctions. So, I am going to keep any junction that 
   has coverage. This will keep the dataset to a manageable size */

data junc_w_cov;
   set junc_counts;
   where apn>0;
   keep event_id;
run;

proc sort data=junc_w_cov nodup;
  by event_id;
proc sort data=junc_counts;
  by event_id;
run;

data junc_counts_dtct;
  merge junc_counts (in=in1) junc_w_cov (in=in2);
  by event_id;
  if in1 and in2;
run;

data event.star_junctions_benchmark_sim;
   set star_juncs;
run;

data event.catalog_junctions_benchmark_sim;
   set junc_counts_dtct;
run;


