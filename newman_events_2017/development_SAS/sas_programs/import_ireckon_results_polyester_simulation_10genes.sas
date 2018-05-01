/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Import iReckon results for the 10 genes with NIC transcripts + MISO-selected genes simulation */

%macro importIR(sample);

proc import datafile="/mnt/store/event_sandbox/polyester_simulated_data/10genes_nic/ireckon_output_10genes_nic_simulation/&sample./result.gtf"
    out=&sample._ireckon dbms=tab replace;
    getnames=no; guessingrows=140000;
run;

/* Parse results GTF */

data &sample._ireckon2;
  set &sample._ireckon;
  length sample_id $10.;
  length gene_id $100.;
  length transcript_id $100.;
  format exon_number best32. ;
  format rpkm best32. ;
  format frac1 best32. ;
  format conf_lo best32. ;
  format conf_hi best32. ;
  format frac2 best32. ;
  format cov best32. ;
  sample_id="&sample.";
  if VAR3="transcript" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
    rpkm=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
    frac1=compress(scan(scan(VAR9,4,";"),2,'"'))+0;
    conf_lo=compress(scan(scan(VAR9,5,";"),2,'"'))+0;
    conf_hi=compress(scan(scan(VAR9,6,";"),2,'"'))+0;
    frac2=compress(scan(scan(VAR9,7,";"),2,'"'))+0;
    cov=compress(scan(scan(VAR9,8,";"),2,'"'))+0;
  end;

  else if VAR3="exon" then do;
    gene_id=compress(scan(scan(VAR9,1,";"),2,'"'));
    transcript_id=compress(scan(scan(VAR9,2,";"),2,'"'));
    exon_number=compress(scan(scan(VAR9,3,";"),2,'"'))+0;
    rpkm=compress(scan(scan(VAR9,4,";"),2,'"'))+0;
    frac1=compress(scan(scan(VAR9,5,";"),2,'"'))+0;
    conf_lo=compress(scan(scan(VAR9,6,";"),2,'"'))+0;
    conf_hi=compress(scan(scan(VAR9,7,";"),2,'"'))+0;
    frac2=compress(scan(scan(VAR9,8,";"),2,'"'))+0;
    cov=compress(scan(scan(VAR9,9,";"),2,'"'))+0;
  end;
  else delete;
  drop VAR9;
  rename VAR1=chr VAR2=source VAR3=feature_type VAR4=start VAR5=stop VAR6=score
         VAR7=strand VAR8=frame;
  run;

%mend;

%importIR(sample_01);
%importIR(sample_02);
%importIR(sample_03);
%importIR(sample_04);
%importIR(sample_05);
%importIR(sample_06);


/* Stack results, subset transcripts and flag novel, IR, etc. */

data all_results;
   set sample_01_ireckon2 sample_02_ireckon2 sample_03_ireckon2
       sample_04_ireckon2 sample_05_ireckon2 sample_06_ireckon2;
   where feature_type="transcript";
  if index(transcript_id,"unspliced") > 0 then flag_novel_unspliced=1; else flaG_novel_unspliced=0;
  if index(transcript_id,"novel") > 0 then flag_novel_isoform=1; else flaG_novel_isoform=0;
  if index(transcript_id,"Intron") > 0 then flag_novel_IR=1; else flaG_novel_IR=0;
  if index(transcript_id,"NM_") > 0 
     or index(transcript_id,"NR_") > 0
     or index(transcript_id,"XM_") > 0
     or index(transcript_id,"XR_") > 0
     then flag_known_isoform=1; else flag_known_isoform=0;
run;


/* Make permenant */

data event.polyester_60gene_sim_ireckon;
  set all_results;
run;



proc freq data= all_results noprint;
  tables sample_id*flaG_known_isoform*flag_novel_isoform*flag_novel_unspliced*flag_novel_IR / out=iso_by_type;
run;

proc print data=iso_by_type;
run;


/*
               flag_      flag_
               known_     novel_    flag_novel_      flag_
 sample_id    isoform    isoform     unspliced     novel_IR    COUNT    PERCENT
 sample_01       0          0            1             0         10      0.3672
 sample_01       0          1            0             0        153      5.6188
 sample_01       0          1            0             1          2      0.0734
 sample_01       1          0            0             0        302     11.0907
 sample_02       0          0            1             0          7      0.2571
 sample_02       0          1            0             0        144      5.2883
 sample_02       0          1            0             1          1      0.0367
 sample_02       1          0            0             0        324     11.8986
 sample_02       1          0            0             1          1      0.0367
 sample_03       0          0            1             0          5      0.1836
 sample_03       0          1            0             0        134      4.9210
 sample_03       0          1            0             1          1      0.0367
 sample_03       1          0            0             0        313     11.4947
 sample_03       1          0            0             1          1      0.0367
 sample_04       0          0            1             0          9      0.3305
 sample_04       0          1            0             0        142      5.2148
 sample_04       0          1            0             1          2      0.0734
 sample_04       1          0            0             0        321     11.7885
 sample_05       0          0            1             0          8      0.2938
 sample_05       0          1            0             0        150      5.5086
 sample_05       1          0            0             0        295     10.8336
 sample_05       1          0            0             1          2      0.0734
 sample_06       0          0            1             0          8      0.2938
 sample_06       0          1            0             0        108      3.9662
 sample_06       0          1            0             1          2      0.0734
 sample_06       1          0            0             0        278     10.2093

(1) Count the number of samples each transcript is observed in
(2) Of the "knowns", flag if from the 10000 transcript list
*/

proc sort data=all_results;
  by gene_id transcript_id;
run;

proc freq data=all_results noprint;
  by gene_id;
  tables transcript_id / out=xs_count;
run;

/* Reflag */

data xs_count_flag;
  set xs_count;
  if index(transcript_id,"unspliced") > 0 then flag_novel_unspliced=1; else flaG_novel_unspliced=0;
  if index(transcript_id,"novel") > 0 then flag_novel_isoform=1; else flaG_novel_isoform=0;
  if index(transcript_id,"Intron") > 0 then flag_novel_IR=1; else flaG_novel_IR=0;
  if index(transcript_id,"NM_") > 0 
     or index(transcript_id,"NR_") > 0
     or index(transcript_id,"XM_") > 0
     or index(transcript_id,"XR_") > 0
     then flag_known_isoform=1; else flag_known_isoform=0;
  if count=6 then flag_all_samples=1; else flaG_all_samples=0;
  if count ge 3 then flag_50perc_samples=1; else flaG_50perc_samples=0;
run;

data sim_list;
  set event.polyester_xs_list_60genes;
run;

proc sort data=xs_count_flag;
  by transcript_id;
proc sort data=sim_list;
  by transcript_id;
run;

data xs_counts_flag_60gn;
  merge xs_count_flag (in=in1) sim_list (in=in2);
  by transcript_id;
  if in2 then flag_simulated_transcript=1; else flag_simulated_transcript=0;
  if in1 then output;
run;

proc freq data=xs_counts_flag_60gn noprint;
  tables flag_known_isoform*flag_novel_isoform*flag_novel_IR*flag_novel_unspliced*
         flag_all_samples*flag_simulated_transcript / out = xs_counts_all_samples;
  tables flag_known_isoform*flag_novel_isoform*flag_novel_IR*flag_novel_unspliced*
         flag_50perc_samples*flag_simulated_transcript / out = xs_counts_50perc_samples;
run;


proc freq data=xs_counts_flag_60gn noprint;
  tables flag_known_isoform*flag_novel_isoform*flag_novel_IR*flag_novel_unspliced*
         flag_all_samples*flag_simulated_transcript / out = xs_counts_all_samples2;
run;
proc freq data=xs_counts_flag_60gn noprint;
  tables flag_known_isoform*flag_novel_isoform*flag_novel_IR*flag_novel_unspliced*
         flag_50perc_samples*flag_simulated_transcript / out = xs_counts_50perc_samples2;
run;

proc print data=xs_counts_all_samples;
run;

proc print data=xs_counts_50perc_samples;
run;

/* ALL:

flag_     flag_
known_    novel_     flag_    flag_novel_   flag_all_   flag_simulated_
soform   isoform   novel_IR    unspliced     samples       transcript     COUNT

  0         0          0           1            0              0            11
  0         0          0           1            1              0             3
  0         1          0           0            0              0           534
  0         1          0           0            1              0             9
  0         1          1           0            0              0             5
  1         0          0           0            0              1           332
  1         0          0           0            1              1           145
  1         0          1           0            0              0             3



50PERC:

  flag_     flag_
  known_    novel_     flag_    flag_novel_   flag_50perc_   flag_simulated_
 isoform   isoform   novel_IR    unspliced       samples        transcript     COUNT

    0         0          0           1              0               0             6
    0         0          0           1              1               0             8
    0         1          0           0              0               0           479
    0         1          0           0              1               0            64
    0         1          1           0              0               0             4
    0         1          1           0              1               0             1
    1         0          0           0              0               1           164
    1         0          0           0              1               1           313
    1         0          1           0              0               0             3


*/


