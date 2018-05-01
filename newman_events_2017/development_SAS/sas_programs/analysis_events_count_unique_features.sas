ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Count detected features by transcript-specificity (unique, common, constitutive) from all features
   I am doing this on the length-filtered list of features
   I already have some of these flags set, but not for splicing events, so I am going to
   remake this */

* First, calculate the number of transcripts per gene;
data xs2gene;
  set refseq.ref_transcript2geneid;
run;

proc freq data=xs2gene noprint;
  tables gene_id / out=xs_per_gene;
run;

data xs_per_gene2;
  set xs_per_gene;
  rename count=xscripts_per_gene;
  drop PERCENT;
run;

/* Merge with feature2xs2gene */
proc sort data=event.feature2xs2gene;
   by gene_id;
proc sort data=xs_per_gene2;
   by gene_id;
run;

data feat2xs2gene;
  merge event.feature2xs2gene (in=in1) xs_per_gene2 (in=in2);
   by gene_id;
  if in1 and in2;
run;

/* Calc transcripts, genes per feature */
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
                  0      925465      100.00        925465       100.00
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

data event.flag_features_by_specificity;
   label genes_per_feature=' ' xscripts_per_feature=' ' max_xscripts_per_feature=' ';
   set flag_features;
run;

/* Counts -- I'll be doing more in python, but I want some counts here too */
proc freq data=event.flag_features_by_specificity;
   tables flag_feature_unique flag_feature_common flag_feature_constitutive
          flag_multigene*flag_feature_common flag_multigene*flag_feature_constitutive;
run;


/*
 flag_feature_                             Cumulative    Cumulative
        unique    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0      650787       70.32        650787        70.32
             1      274678       29.68        925465       100.00


 flag_feature_                             Cumulative    Cumulative
        common    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0      670830       72.49        670830        72.49
             1      254635       27.51        925465       100.00


 flag_feature_                             Cumulative    Cumulative
  constitutive    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0      529313       57.19        529313        57.19
             1      396152       42.81        925465       100.00


 Table of flag_multigene by flag_feature_common

      flag_multigene
                flag_feature_common

      Frequency|
      Percent  |
      Row Pct  |
      Col Pct  |       0|       1|  Total
      ---------+--------+--------+
             0 | 661475 | 241917 | 903392
               |  71.47 |  26.14 |  97.61
               |  73.22 |  26.78 |
               |  98.61 |  95.01 |
      ---------+--------+--------+
             1 |   9355 |  12718 |  22073
               |   1.01 |   1.37 |   2.39
               |  42.38 |  57.62 |
               |   1.39 |   4.99 |
      ---------+--------+--------+
      Total      670830   254635   925465
                  72.49    27.51   100.00

 Table of flag_multigene by flag_feature_constitutive

         flag_multigene
                   flag_feature_constitutive

         Frequency|
         Percent  |
         Row Pct  |
         Col Pct  |       0|       1|  Total
         ---------+--------+--------+
                0 | 516595 | 386797 | 903392
                  |  55.82 |  41.79 |  97.61
                  |  57.18 |  42.82 |
                  |  97.60 |  97.64 |
         ---------+--------+--------+
                1 |  12718 |   9355 |  22073
                  |   1.37 |   1.01 |   2.39
                  |  57.62 |  42.38 |
                  |   2.40 |   2.36 |
         ---------+--------+--------+
         Total      529313   396152   925465
                     57.19    42.81   100.00
*/


