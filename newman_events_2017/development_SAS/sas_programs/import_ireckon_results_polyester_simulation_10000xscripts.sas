/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Import iReckon results for 10000 transcript simulation */

%macro importIR(sample);

proc import datafile="!MCLAB/event_analysis/alignment_output/ireckon/polyester_10000_xscripts/&sample..gtf"
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

data event.polyester_10k_simulation_ireckon;
  set all_results;
run;


proc freq data= all_results noprint;
  tables sample_id*flaG_known_isoform*flag_novel_isoform*flag_novel_unspliced*flag_novel_IR / out=iso_by_type;
run;

proc print data=iso_by_type;
run;

/* For sample 1:
               flag_      flag_
               known_     novel_    flag_novel_      flag_
 sample_id    isoform    isoform     unspliced     novel_IR    COUNT    PERCENT

 sample_01       0          0            1             0         481     0.6524
 sample_01       0          1            0             0        2147     2.9119
 sample_01       0          1            0             1          13     0.0176
 sample_01       1          0            0             0        9569    12.9781
 sample_01       1          0            0             1          13     0.0176


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

data xs10000;
  set event.polyester_xs_list_10k;
run;

proc sort data=xs_count_flag;
  by transcript_id;
proc sort data=xs10000;
  by transcript_id;
run;

data xs_counts_flag_10k;
  merge xs_count_flag (in=in1) xs10000 (in=in2);
  by transcript_id;
  if in2 then flag_simulated_transcript=1; else flag_simulated_transcript=0;
  if in1 then output;
run;

proc freq data=xs_counts_flag_10k noprint;
  tables flag_known_isoform*flag_novel_isoform*flag_novel_IR*flag_novel_unspliced*
         flag_all_samples*flag_simulated_transcript / out = xs_counts_all_samples;
  tables flag_known_isoform*flag_novel_isoform*flag_novel_IR*flag_novel_unspliced*
         flag_50perc_samples*flag_simulated_transcript / out = xs_counts_50perc_samples;
run;

proc print data=xs_counts_all_samples;
run;

proc print data=xs_counts_50perc_samples;
run;

/*
 flag_     flag_
 known_    novel_     flag_    flag_novel_   flag_all_   flag_simulated_
isoform   isoform   novel_IR    unspliced     samples       transcript     COUNT

   0         0          0           1            0              0            792
   0         0          0           1            1              0            224
   0         1          0           0            0              0           5393
   0         1          0           0            1              0            499
   0         1          1           0            0              0             48
   0         1          1           0            1              0              2
   1         0          0           0            0              0           1565
   1         0          0           0            0              1           2141
   1         0          0           0            1              0           1197
   1         0          0           0            1              1           6370
   1         0          1           0            0              0             29
   1         0          1           0            1              0              2


 flag_     flag_
 known_    novel_     flag_    flag_novel_   flag_50perc_   flag_simulated_
isoform   isoform   novel_IR    unspliced       samples        transcript     COUNT

   0         0          0           1              0               0            581
   0         0          0           1              1               0            435
   0         1          0           0              0               0           4014
   0         1          0           0              1               0           1878
   0         1          1           0              0               0             41
   0         1          1           0              1               0              9
   1         0          0           0              0               0            917
   1         0          0           0              0               1            731
   1         0          0           0              1               0           1845
   1         0          0           0              1               1           7780
   1         0          1           0              0               0             22
   1         0          1           0              1               0              9
*/


