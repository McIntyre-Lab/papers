
ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Comparing the updated feature specificity flags with the flags pre-gene-filtering to examine any shifts in transcript specificity of features. This will be on all features together (fusions, fragments, splicing), as I am only interested in checking if there are substantial overall changes */


data old_flags;
  set event.flag_features_by_specificity;
  keep feature_id flag_multigene flag_feature_unique flag_feature_common flag_feature_constitutive;
  rename flag_multigene=flag_multigene_old
         flag_feature_unique=flag_unique_old
         flag_feature_common=flag_common_old
         flag_feature_constitutive=flag_constitutive_old;
run;


data new_flags;
  set event.flag_features_by_specificity_exp;
  keep feature_id flag_multigene flag_feature_unique flag_feature_common flag_feature_constitutive;
  rename flag_multigene=flag_multigene_new
         flag_feature_unique=flag_unique_new
         flag_feature_common=flag_common_new
         flag_feature_constitutive=flag_constitutive_new;
run;

proc sort data=old_flags;
   by feature_id;
proc sort data=new_flags;
   by feature_id;
run;

data old_and_new_flags;
  merge old_flags (in=in1) new_flags (in=in2);
  by feature_id;
  if in1 and in2;
run;

/* Count: how many multigene features became single gene? */
proc freq data=old_and_new_flags;
   tables flag_multigene_old*flag_multigene_new;
run;

/*

  flag_multigene_old
            flag_multigene_new

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 | 793713 |      0 | 793713
           |  98.13 |   0.00 |  98.13
           | 100.00 |   0.00 |
           |  99.04 |   0.00 |
  ---------+--------+--------+
         1 |   7665 |   7482 |  15147
           |   0.95 |   0.93 |   1.87
           |  50.60 |  49.40 |
           |   0.96 | 100.00 |
  ---------+--------+--------+
  Total      801378     7482   808860
              99.07     0.93   100.00


*/

/* Count: compare new unique to old flags -- expect "old unique" to remain unique */
proc freq data=old_and_new_flags noprint;
   tables flag_unique_new*flag_unique_old*flag_common_old*flag_constitutive_old / out=uniq_check;
run;

proc print data=uniq_check;
run;

/*
  flag_      flag_      flag_
unique_    unique_    common_    flag_constitutive_
  new        old        old              old            COUNT    PERCENT

   0          0          0                1            349502    43.2092
   0          0          1                0            227422    28.1164
   1          0          0                1              3348     0.4139
   1          0          1                0              3405     0.4210
   1          1          0                0            225183    27.8396

*/



/* Count: compare new common to old flags */
proc freq data=old_and_new_flags noprint;
   tables flag_common_new*flag_unique_old*flag_common_old*flag_constitutive_old / out=common_check;
run;

proc print data=common_check;
run;

/*
  flag_      flag_      flag_
 common_    unique_    common_    flag_constitutive_
   new        old        old              old            COUNT    PERCENT

    0          0          0                1            352850    43.6231
    0          0          1                0              8365     1.0342
    0          1          0                0            225183    27.8396
    1          0          1                0            222462    27.5032

*/


/* Count: compare new constitutive to old flags */
proc freq data=old_and_new_flags noprint;
   tables flag_constitutive_new*flag_unique_old*flag_common_old*flag_constitutive_old / out=constit_check;
run;

proc print data=constit_check;
run;


/*
                        flag_      flag_
 flag_constitutive_    unique_    common_    flag_constitutive_
         new             old        old              old            COUNT    PERCENT

          0               0          0                1              3348     0.4139
          0               0          1                0            225867    27.9241
          0               1          0                0            225183    27.8396
          1               0          0                1            349502    43.2092
          1               0          1                0              4960     0.6132

*/
