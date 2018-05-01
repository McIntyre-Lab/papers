ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Using the update feature2transcript list, update feature specificity to consider only those transcripts from genes that are expressed. This shouldn't change what is unique/common/constitutive
too much, and we can check the new classifications against the old to see what has shifted

E.g. If a fragment initially was assigned to two different genes and one of those genes is not expressed, then fragment is now only assigned to expressed gene */

* First, calculate the number of transcripts per gene;

data xs2gene;
   set event.feature2xs2gene_exp_only;
   keep gene_id transcript_id;
run;

proc sort data=xs2gene nodup;
  by gene_id transcript_id;
proc freq data=xs2gene noprint;
   tables gene_id / out=xs_per_gene;
run;

data xs_per_gene2;
  set xs_per_gene;
  rename count=xscripts_per_gene;
  drop PERCENT;
run;

/* Merge with feature2xs2gene */
proc sort data=event.feature2xs2gene_exp_only nodup;
   by gene_id;
proc sort data=xs_per_gene2;
   by gene_id;
run;

data feat2xs2gene;
  merge event.feature2xs2gene (in=in1) xs_per_gene2 (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* Calc transcripts and genes per features */

data feat2xs;
  set feat2xs2gene;
  keep feature_id transcript_id ;
run;

proc sort data=feat2xs nodup;
  by feature_id transcript_id;
proc freq data=feat2xs noprint;
  tables feature_id / out=xs_per_feature;
run;

data xs_per_feature2;
  set xs_per_feature;
  drop percent;
  rename count=xscripts_per_feature;
run;

data feat2gene;
  set feat2xs2gene;
  drop transcript_id;
run;

proc sort data=feat2gene nodup;
   by feature_id gene_id xscripts_per_gene;
proc means data=feat2gene noprint;
   by feature_id;
   var xscripts_per_gene;
   output out=max_xscripts_per_feat sum=max_xscripts_per_feature;
proc freq data=feat2gene noprint;
   tables feature_id / out=genes_per_feature;
run;

data genes_per_feature2;
  set genes_per_feature;
  drop percent;
  rename count=genes_per_feature;
run;

data max_xscripts_per_feat2;
  set max_xscripts_per_feat;
  keep feature_id max_xscripts_per_feature;
run;


/* Merge all together, make flags */

proc sort data=genes_per_feature2;
  by feature_id;
proc sort data=xs_per_feature2;
  by feature_id;
proc sort data=max_xscripts_per_feat2;
  by feature_id;
run;

data feat_w_gene_xs_counts;
  merge genes_per_feature2 xs_per_feature2 max_xscripts_per_feat2;
  by feature_id;
run;

/* Check to make sure that the number of xscripts per feature does not exceed the
   maximum possible of xscripts*/

data check_xs_num;
  set feat_w_gene_xs_counts;
  if xscripts_per_feature > max_xscripts_per_feature then flag_feat2xs_bad=1;
  else flag_feat2xs_bad=0;
run;

proc freq data=check_xs_num;
  tables flag_feat2xs_bad;
run;

/*
                                               Cumulative    Cumulative
  flag_feat2xs_bad    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0      814627      100.00        814627       100.00

*/

/*check 2: are there features with multiple genes, but one transcript? there should be 0 */

data check;
  set feat_w_gene_xs_counts;
  where xscripts_per_feature=1 and genes_per_feature > 1;
run; *0 observations, good!;

/* Set flags for unique, common, constitutive */

data flag_features;
  set feat_W_gene_xs_counts;
  if genes_per_feature = 1 then do;
     flag_multigene=0;
     if xscripts_per_feature=1 then do;
        flag_feature_unique=1; flag_feature_common=0; flag_feature_constitutive=0; end;
     else if xscripts_per_feature=max_xscripts_per_feature then do;
        flag_feature_unique=0; flag_feature_common=0; flag_feature_constitutive=1; end;
     else if xscripts_per_feature < max_xscripts_per_feature then do;
        flag_feature_unique=0; flag_feature_common=1; flag_feature_constitutive=0; end;
     else flag_oops=1;
     end;
   else do;
     flag_multigene=1;
     if xscripts_per_feature=max_xscripts_per_feature then do;
        flag_feature_unique=0; flag_feature_common=0; flag_feature_constitutive=1; end;
     else if xscripts_per_feature < max_xscripts_per_feature then do;
        flag_feature_unique=0; flag_feature_common=1; flag_feature_constitutive=0; end;
     else flag_oops=1;
     end;
run;

/* Make permenant */

data event.flag_features_by_specificity_exp;
   label genes_per_feature=' ' xscripts_per_feature=' ' max_xscripts_per_feature=' ';
   set flag_features;
run;

/* Counts -- I'll be doing more in python, but I want some counts here too */
proc freq data=event.flag_features_by_specificity_exp;
   tables flag_feature_unique flag_feature_common flag_feature_constitutive
          flag_multigene*flag_feature_common flag_multigene*flag_feature_constitutive;
run;


/*


      flag_feature_                             Cumulative    Cumulative
             unique    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0      584271       71.72        584271        71.72
                  1      230356       28.28        814627       100.00


      flag_feature_                             Cumulative    Cumulative
             common    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0      582481       71.50        582481        71.50
                  1      232146       28.50        814627       100.00


      flag_feature_                             Cumulative    Cumulative
       constitutive    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0      462502       56.77        462502        56.77
                  1      352125       43.23        814627       100.00


   flag_multigene
             flag_feature_common

   Frequency|
   Percent  |
   Row Pct  |
   Col Pct  |       0|       1|  Total
   ---------+--------+--------+
          0 | 580521 | 226503 | 807024
            |  71.26 |  27.80 |  99.07
            |  71.93 |  28.07 |
            |  99.66 |  97.57 |
   ---------+--------+--------+
          1 |   1960 |   5643 |   7603
            |   0.24 |   0.69 |   0.93
            |  25.78 |  74.22 |
            |   0.34 |   2.43 |
   ---------+--------+--------+
   Total      582481   232146   814627
               71.50    28.50   100.00

    flag_multigene
              flag_feature_constitutive

    Frequency|
    Percent  |
    Row Pct  |
    Col Pct  |       0|       1|  Total
    ---------+--------+--------+
           0 | 456859 | 350165 | 807024
             |  56.08 |  42.98 |  99.07
             |  56.61 |  43.39 |
             |  98.78 |  99.44 |
    ---------+--------+--------+
           1 |   5643 |   1960 |   7603
             |   0.69 |   0.24 |   0.93
             |  74.22 |  25.78 |
             |   1.22 |   0.56 |
    ---------+--------+--------+
    Total      462502   352125   814627
                56.77    43.23   100.00




*/


