ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Prep and export the following data:
   Dataset (1) : perc_unique_dtct, perc_total_dtct (plus transcript "bins")
		 Make these stacked barplots!!!
   
   Dataset (2) : num_xscripts_per_gene for each "bin" of perc total dtct */

/* Dataset (1) : e.g. 100% : no uniq, 0%, 0-25%, 25-50%, 50-75%, 75-100%, 100% */

data bin_xscripts;
   set event.xscripts_w_unique_by_bin;
   length bin_xscript_perc_total_dtct $8.;
   if perc_features_dtct=1 then bin_xscript_perc_total_dtct="100%";
   else if perc_features_dtct ge 0.75 then bin_xscript_perc_total_dtct="75-100%";
   else if perc_features_dtct ge 0.5 then bin_xscript_perc_total_dtct="50-75%";
   else if perc_features_dtct ge 0.25 then bin_xscript_perc_total_dtct="25-50%";
   else if perc_features_dtct gt 0 then bin_xscript_perc_total_dtct="<25%";
   else delete;
   keep transcript_id 
        bin_xscript_perc_total_dtct bin_xscript_perc_uniq_dtct 
        perc_unique_features_dtct perc_features_dtct;
run;

proc freq data=bin_xscripts;
  tables bin_xscript_perc_total_dtct*bin_xscript_perc_uniq_dtct ;
run;


/*
bin_xscript_perc_total_dtct     bin_xscript_perc_uniq_dtct

Frequency|
Percent  |
Row Pct  |
Col Pct  |0%      |0-25%   |100%    |25-50%  |50-75%  |75-100% |no uniqu|  Total
         |        |        |        |        |        |        |e       |
---------+--------+--------+--------+--------+--------+--------+--------+
100%     |      0 |      0 |   9730 |      0 |      0 |      0 |   5004 |  14734
         |   0.00 |   0.00 |  13.23 |   0.00 |   0.00 |   0.00 |   6.80 |  20.04
         |   0.00 |   0.00 |  66.04 |   0.00 |   0.00 |   0.00 |  33.96 |
         |   0.00 |   0.00 |  76.64 |   0.00 |   0.00 |   0.00 |  14.34 |
---------+--------+--------+--------+--------+--------+--------+--------+
75-100%  |   5117 |     10 |   2240 |    455 |   3331 |    646 |   8089 |  19888
         |   6.96 |   0.01 |   3.05 |   0.62 |   4.53 |   0.88 |  11.00 |  27.05
         |  25.73 |   0.05 |  11.26 |   2.29 |  16.75 |   3.25 |  40.67 |
         |  31.36 |   0.77 |  17.64 |  22.25 |  59.42 |  96.56 |  23.18 |
---------+--------+--------+--------+--------+--------+--------+--------+
50-75%   |   2998 |     16 |    454 |    293 |   1643 |     20 |   5837 |  11261
         |   4.08 |   0.02 |   0.62 |   0.40 |   2.23 |   0.03 |   7.94 |  15.31
         |  26.62 |   0.14 |   4.03 |   2.60 |  14.59 |   0.18 |  51.83 |
         |  18.37 |   1.23 |   3.58 |  14.33 |  29.31 |   2.99 |  16.72 |
---------+--------+--------+--------+--------+--------+--------+--------+
25-50%   |   3207 |     32 |    213 |   1180 |    428 |      1 |   6628 |  11689
         |   4.36 |   0.04 |   0.29 |   1.60 |   0.58 |   0.00 |   9.01 |  15.90
         |  27.44 |   0.27 |   1.82 |  10.09 |   3.66 |   0.01 |  56.70 |
         |  19.66 |   2.46 |   1.68 |  57.70 |   7.63 |   0.15 |  18.99 |
---------+--------+--------+--------+--------+--------+--------+--------+
<25%     |   4994 |   1245 |     59 |    117 |    204 |      2 |   9342 |  15963
         |   6.79 |   1.69 |   0.08 |   0.16 |   0.28 |   0.00 |  12.70 |  21.71
         |  31.28 |   7.80 |   0.37 |   0.73 |   1.28 |   0.01 |  58.52 |
         |  30.61 |  95.55 |   0.46 |   5.72 |   3.64 |   0.30 |  26.77 |
---------+--------+--------+--------+--------+--------+--------+--------+
Total       16316     1303    12696     2045     5606      669    34900    73535
            22.19     1.77    17.27     2.78     7.62     0.91    47.46   100.00

*/


/* For each gene, count the number of transcripts at each "total detected" threshold */


data xs2gene;
  set event.feature2xs2gene_exp_only_nomulti;
  keep transcript_id gene_id;
run;

proc sort data=xs2gene nodup;
   by transcript_id gene_id;
run;

%macro xsPerGene(minDtct,percDtct);

data xscripts_kept;
  set bin_xscripts;
  where perc_features_dtct ge &minDtct.;
  keep transcript_id;
run;

proc sort data=xscripts_kept;
   by transcript_id;
run;

data xs2gene_kept;
  merge xs2gene (in=in1) xscripts_kept (in=in2);
  by transcript_id;
  if in1 and in2;
run;

proc sort data=xs2gene_kept;
   by gene_id;
proc freq data=xs2gene_kept noprint;
   tables gene_id / out=xs_per_gene;
run;

data xs_per_gene_&percDtct.;
   set xs_per_gene;
   keep gene_id count;
   rename count=num_xscript_min_perc_on_&percDtct.;
run;

%mend;

%xsPerGene(1,100);
%xsPerGene(0.75,75);
%xsPerGene(0.5,50);
%xsPerGene(0.25,25);
%xsPerGene(0,0);

proc sort data=xs_per_gene_100;
   by gene_id;
proc sort data=xs_per_gene_75;
   by gene_id;
proc sort data=xs_per_gene_50;
   by gene_id;
proc sort data=xs_per_gene_25;
   by gene_id;
proc sort data=xs_per_gene_0;
   by gene_id;
run;

data xs_per_gene_all;
  merge xs_per_gene_100 (in=in1) xs_per_gene_75 (in=in2) xs_per_gene_50 (in=in3) 
        xs_per_gene_25 (in=in4) xs_per_gene_0 (in=in5) ;
  by gene_id;
  if not in1 then num_xscript_min_perc_on_100=0;
  if not in2 then num_xscript_min_perc_on_75=0;
  if not in3 then num_xscript_min_perc_on_50=0;
  if not in4 then num_xscript_min_perc_on_25=0;
  if not in5 then num_xscript_min_perc_on_0=0;
run;

/* Make permenant */

data event.xscripts_w_unique_by_bin_total;
    set bin_xscripts;
run;


data event.xscripts_per_gene_by_total_dtct;
   set xs_per_gene_all;
run;


/* Export data for plots */

proc export data=event.xscripts_per_gene_by_total_dtct
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_xscripts_per_gene_perc_total_detected_nomulti.csv"
     dbms=csv replace;
run;

proc export data=event.xscripts_w_unique_by_bin_total
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_xscript_bins_uniq_total_detected_nomulti.csv"
     dbms=csv replace;
run;


