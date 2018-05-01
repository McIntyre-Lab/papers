/* Libraries */
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* Compare number of reads simulated per transcript against the proportion of events detected
   There are ~373 transcripts with no coverage that I want to figure out why

   I expect that these transcripts have very few reads on average
*/

/* Import simulated read counts per sample */

%macro importCount(sample);

proc import datafile="/mnt/store/event_sandbox/polyester_simulated_data/10000_refseq_transcripts/&sample._read_counts.txt"
   out=&sample._reads dbms=tab replace;
   getnames=no; guessingrows=10000;
run;

data &sample._reads2;
  length sample_id $9.;
  format read_count best32.;
  length transcript_id $20.;
  set &sample._reads;
  sample_id="&sample.";
  read_count=scan(VAR1,1," ") + 0;
  transcript_id=compress(scan(VAR1,2," "));
  drop VAR1;
run;

%mend;

%importCount(sample_01);
%importCount(sample_02);
%importCount(sample_03);
%importCount(sample_04);
%importCount(sample_05);
%importCount(sample_06);


data read_counts_by_sample;
   set sample_01_reads2 sample_02_reads2 sample_03_reads2
       sample_04_reads2 sample_05_reads2 sample_06_reads2;
   log_read_count=log(read_count+1);
run;

proc sort data=read_counts_by_sample;
  by transcript_id;
proc means data=read_counts_by_sample noprint;
  by transcript_id;
  var log_read_count;
  output out=mean_read_count_by_xs min=min_count max=max_count mean=mean_count;
run;

data perc_dtct;
   set event.bin_xs_by_dtct_apn0_10k;
   keep perc_features_dtct perc_junctions_dtct perc_fragments_dtct transcript_id;
run;

proc sort data=mean_read_count_by_xs;
  by transcript_id;
proc sort data=perc_dtct;
  by transcript_id;
run;

data read_count2perc_dtct;
  merge mean_read_count_by_xs (in=in1) perc_dtct (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc corr data=read_count2perc_dtct pearson;
   var min_count max_count mean_count perc_features_dtct perc_junctions_dtct perc_fragments_dtct;
run;

/*

                                                                 perc_         perc_         perc_
                                                    mean_    features_    junctions_    fragments_
                       min_count    max_count       count         dtct          dtct          dtct

min_count                1.00000      0.99897     0.99965      0.17027       0.16808       0.12647
                                       <.0001      <.0001       <.0001        <.0001        <.0001
                           10000        10000       10000        10000          9482         10000

max_count                0.99897      1.00000     0.99968      0.16875       0.16361       0.12554
                          <.0001                   <.0001       <.0001        <.0001        <.0001
                           10000        10000       10000        10000          9482         10000

mean_count               0.99965      0.99968     1.00000      0.16979       0.16614       0.12621
                          <.0001       <.0001                   <.0001        <.0001        <.0001
                           10000        10000       10000        10000          9482         10000

perc_junctions_dtct      0.16808      0.16361     0.16614      0.37029       1.00000       0.17779
                          <.0001       <.0001      <.0001       <.0001                      <.0001
                            9482         9482        9482         9482          9482          9482

perc_fragments_dtct      0.12647      0.12554     0.12621      0.97327       0.17779       1.00000
                          <.0001       <.0001      <.0001       <.0001        <.0001
                           10000        10000       10000        10000          9482         10000


Log read count is poorly correlated with the proportion of features detected

Plot mean read count against perc_features_dtct
*/

proc sgplot data=read_count2perc_dtct;
   scatter x=mean_count  y=perc_features_dtct;
run;

