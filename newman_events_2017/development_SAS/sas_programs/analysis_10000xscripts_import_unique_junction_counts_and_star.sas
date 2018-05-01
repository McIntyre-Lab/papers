/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* For the 10000-transcript simulation, import the STAR junction results and EA junction coverage counts and check
   to see which method captures the junctions of these genes */

/* Import STAR */

%macro importSTAR(sample);

proc import datafile="!MCLAB/event_analysis/alignment_output/aln_mm10_star_10k_simulation/&sample.SJ.out.tab"
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

%importSTAR(sample_01);
%importSTAR(sample_02);
%importSTAR(sample_03);
%importSTAR(sample_04);
%importSTAR(sample_05);
%importSTAR(sample_06);

/* Import coverage counts */

%macro importCC(sample);

    data &sample._counts    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile "!MCLAB/event_analysis/alignment_output/coverage_counts_splicing_10k/cvrg_cnts_&sample..csv" delimiter
= ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat sample_id $9. ;
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
       format sample_id $9. ;
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

%importCC(sample_01);
%importCC(sample_02);
%importCC(sample_03);
%importCC(sample_04);
%importCC(sample_05);
%importCC(sample_06);

/* Stack outputs and save */

data star_juncs;
   set sample_01_star2 sample_02_star2 sample_03_star2
       sample_04_star2 sample_05_star2 sample_06_star2;
run;

data junc_counts;
   set sample_01_counts sample_02_counts sample_03_counts
       sample_04_counts sample_05_counts sample_06_counts;
run;

/* For catalog junctions, I really only want to keep the junctions from the transcripts selected, but I also want
   to double check that reads didn't got to the wrong junctions. So, I am going to keep any junction that 
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

data event.star_junctions_10000xs;
   set star_juncs;
run;

data event.catalog_junctions_10000xs;
   set junc_counts_dtct;
run;


