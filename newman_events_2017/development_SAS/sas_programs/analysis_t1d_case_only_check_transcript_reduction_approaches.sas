data reduce_1;
  set eventloc.xs_reduced_by_tpm_subj_cell_v1;
  where flag_tpm_reduced=1;
  keep transcript_id;
run;


data reduce_2;
  set eventloc.xs_reduced_by_tpm_subj_cell_v2;
  where flag_tpm_reduced=1;
  keep transcript_id;
run;


proc sort data=reduce_1 nodup;
  by transcript_id;
run;

proc sort data=reduce_2 nodup;
  by transcript_id;
run;

data check_reduce;
  merge reduce_1 (in=in1) reduce_2 (in=in2);
  by transcript_id;
  if in1 then flag_reduced_by_subj=1; else flag_reduced_by_subj=0;
  if in2 then flag_reduced_by_mean=1; else flag_reduced_by_mean=0;
run;


proc freq data=check_reduce;
   table flag_reduced_by_subj*flag_reduced_by_mean;
run;

/*

  flag_reduced_by_subj
            flag_reduced_by_mean

  Frequency|
  Percent  |
  Row Pct  |
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         1 |   3392 |  22401 |  25793
           |  13.15 |  86.85 | 100.00
           |  13.15 |  86.85 |
           | 100.00 | 100.00 |
  ---------+--------+--------+
  Total        3392    22401    25793
              13.15    86.85   100.00


All isofomrs in the reduced set based off the mean TPM are in the set that are in the reduced set based off the TPM by subject. So now I want to see of the ones only in the "by-subject", how many have counts in all subjects? (these should not have all ~241 samples..)
*/




