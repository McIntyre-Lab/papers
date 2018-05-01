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
Frequency|
Percent  |
Row Pct  |
Col Pct  |0%      |0-25%   |100%    |25-50%  |50-75%  |75-100% |no uniqu|  Total
         |        |        |        |        |        |        |e       |
---------+--------+--------+--------+--------+--------+--------+--------+
100%     |      0 |      0 |  12421 |      0 |      0 |      0 |   6917 |  19338
         |   0.00 |   0.00 |  12.21 |   0.00 |   0.00 |   0.00 |   6.80 |  19.01
         |   0.00 |   0.00 |  64.23 |   0.00 |   0.00 |   0.00 |  35.77 |
         |   0.00 |   0.00 |  72.85 |   0.00 |   0.00 |   0.00 |  14.47 |
---------+--------+--------+--------+--------+--------+--------+--------+
75-100%  |   7432 |     15 |   3472 |    650 |   4802 |    939 |  12265 |  29575
         |   7.31 |   0.01 |   3.41 |   0.64 |   4.72 |   0.92 |  12.06 |  29.07
         |  25.13 |   0.05 |  11.74 |   2.20 |  16.24 |   3.17 |  41.47 |
         |  31.86 |   0.90 |  20.36 |  23.44 |  59.21 |  95.33 |  25.66 |
---------+--------+--------+--------+--------+--------+--------+--------+
50-75%   |   4568 |     38 |    750 |    501 |   2356 |     41 |   8245 |  16499
         |   4.49 |   0.04 |   0.74 |   0.49 |   2.32 |   0.04 |   8.11 |  16.22
         |  27.69 |   0.23 |   4.55 |   3.04 |  14.28 |   0.25 |  49.97 |
         |  19.58 |   2.28 |   4.40 |  18.07 |  29.05 |   4.16 |  17.25 |
---------+--------+--------+--------+--------+--------+--------+--------+
25-50%   |   4536 |     82 |    314 |   1418 |    643 |      3 |   8757 |  15753
         |   4.46 |   0.08 |   0.31 |   1.39 |   0.63 |   0.00 |   8.61 |  15.49
         |  28.79 |   0.52 |   1.99 |   9.00 |   4.08 |   0.02 |  55.59 |
         |  19.44 |   4.92 |   1.84 |  51.14 |   7.93 |   0.30 |  18.32 |
---------+--------+--------+--------+--------+--------+--------+--------+
<25%     |   6793 |   1533 |     93 |    204 |    309 |      2 |  11621 |  20555
         |   6.68 |   1.51 |   0.09 |   0.20 |   0.30 |   0.00 |  11.42 |  20.21
         |  33.05 |   7.46 |   0.45 |   0.99 |   1.50 |   0.01 |  56.54 |
         |  29.12 |  91.91 |   0.55 |   7.36 |   3.81 |   0.20 |  24.31 |
---------+--------+--------+--------+--------+--------+--------+--------+
Total       23329     1668    17050     2773     8110      985    47805   101720
            22.93     1.64    16.76     2.73     7.97     0.97    47.00   100.00

*/


/* For each gene, count the number of transcripts at each "total detected" threshold */


data xs2gene;
  set event.feature2xs2gene;
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
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_xscripts_per_gene_perc_total_detected.csv"
     dbms=csv replace;
run;

proc export data=event.xscripts_w_unique_by_bin_total
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_xscript_bins_uniq_total_detected.csv"
     dbms=csv replace;
run;


