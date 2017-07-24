
ods listing; ods html close;

libname event '!MCLAB/event_analysis/sas_data';

/* get APN<0/5 50% and 75% lists */

data list1 list2;
  set event.xscripts_w_unique_by_bin;
  if perc_features_dtct >= 0.5 then output list1;
  if perc_features_dtct >= 0.75 then output list2;
   keep transcript_id;
run;

data list3 list4;
  set event.bin_xscripts_by_dtct_apn5;
  if perc_features_dtct >= 0.5 then output list3;
  if perc_features_dtct >= 0.75 then output list4;
   keep transcript_id;
run;




data xs2gene;
  set event.feature2xs2gene_exp_only_nomulti;
  keep transcript_id gene_id;
  run;

proc sort data=xs2gene nodup;
  by transcript_id gene_id;
proc sort data=list1;
  by transcript_id;
proc sort data=list2;
  by transcript_id;
proc sort data=list3;
  by transcript_id;
proc sort data=list4;
  by transcript_id;
  run;

  data list1_w_gene;
     merge xs2gene (in=in1) list1 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;

  data list2_w_gene;
     merge xs2gene (in=in1) list2 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;


  data list3_w_gene;
     merge xs2gene (in=in1) list3 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;

  data list4_w_gene;
     merge xs2gene (in=in1) list4 (in=in2);
     by transcript_id;
     if in1 and in2;
     run;




    proc export data=list1_w_gene
       outfile="!MCLAB/event_analysis/pacbio_plus_apn0_50perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;

    proc export data=list2_w_gene
       outfile="!MCLAB/event_analysis/pacbio_plus_apn0_75perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;
    proc export data=list3_w_gene
       outfile="!MCLAB/event_analysis/pacbio_plus_apn5_50perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;
    proc export data=list4_w_gene
       outfile="!MCLAB/event_analysis/pacbio_plus_apn5_75perc_gene2xs.txt"
       dbms=tab replace;
       putnames=no;
       run;

