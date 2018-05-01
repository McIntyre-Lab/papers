ods listing; ods html close;
libname event '!MCLAB/event_analysis/sas_data';

data bin_xscripts;
   set event.hg19_flag_xscripts_w_unique;
   length bin_xscript_perc_total_dtct $8.;
   max_perc_features_dtct=max(perc_features_dtct_cd4,perc_features_dtct_cd8,perc_features_dtct_cd19);
   if max_perc_features_dtct=1 then bin_xscript_perc_total_dtct="100%";
   else if max_perc_features_dtct ge 0.75 then bin_xscript_perc_total_dtct="75-100%";
   else if max_perc_features_dtct ge 0.5 then bin_xscript_perc_total_dtct="50-75%";
   else if max_perc_features_dtct ge 0.25 then bin_xscript_perc_total_dtct="25-50%";
   else if max_perc_features_dtct gt 0 then bin_xscript_perc_total_dtct="<25%";
   else delete;
   keep transcript_id 
        bin_xscript_perc_total_dtct max_perc_features_dtct;
run;


/* For each gene, count the number of transcripts at each "total detected" threshold */


data xs2gene;
  set event.hg19_feature2xs2gene_exp_only2;
  keep transcript_id gene_id;
run;

proc sort data=xs2gene nodup;
   by transcript_id gene_id;
run;

%macro xsPerGene(minDtct,percDtct);

data xscripts_kept;
  set bin_xscripts;
  where max_perc_features_dtct ge &minDtct.;
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

data event.t1d_xs_per_gene_by_total_dtct;
   set xs_per_gene_all;
run;


/* Export data for plots */

proc export data=event.t1d_xs_per_gene_by_total_dtct
     outfile="!MCLAB/event_analysis/analysis_output/event_analysis_t1d_xscripts_per_gene_perc_total_detected_nomulti.csv"
     dbms=csv replace;
run;


