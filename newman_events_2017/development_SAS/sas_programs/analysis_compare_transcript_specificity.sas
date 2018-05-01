
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
         0 | 799222 |      0 | 799222
           |  98.11 |   0.00 |  98.11
           | 100.00 |   0.00 |
           |  99.03 |   0.00 |
  ---------+--------+--------+
         1 |   7802 |   7603 |  15405
           |   0.96 |   0.93 |   1.89
           |  50.65 |  49.35 |
           |   0.97 | 100.00 |
  ---------+--------+--------+
  Total      807024     7603   814627
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

   0          0          0                1            351811    43.1868
   0          0          1                0            232460    28.5358
   1          0          0                1              1039     0.1275
   1          0          1                0              1665     0.2044
   1          1          0                0            227652    27.9456



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

    0          0          0                1            352850    43.3143
    0          0          1                0              1979     0.2429
    0          1          0                0            227652    27.9456
    1          0          1                0            232146    28.4972


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

         0               0          0                1              1039     0.1275
         0               0          1                0            233811    28.7016
         0               1          0                0            227652    27.9456
         1               0          0                1            351811    43.1868
         1               0          1                0               314     0.0385

*/
