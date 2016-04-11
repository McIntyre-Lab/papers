/* Calculate means for each treatment by fusion?*/

libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data" ;

data prep_means;
  set ribo.all_coverage_counts_with_key;
  where trt ne "";
  run;

proc sort data=prep_means;
  by fusion_id trt;
  run;

/*Calculate mean and SD for logRPKM and transpose*/ 
proc means data=prep_means noprint;
  var logrpkm;
  by fusion_id;
  class trt;
  output out=logrpkm_mean mean=mean stddev=SD;
  run;

data logrpkm_mean;
  set logrpkm_mean;
  where _type_ = 1;
  run;

proc transpose data=logrpkm_mean out=flip_logrpkm_mean suffix=_mean_logrpkm;
  by fusion_id;
  var mean;
  id trt;
  run;

proc transpose data=logrpkm_mean out=flip_logrpkm_sd suffix=_SD_logrpkm;
  by fusion_id;
  var SD;
  id trt;
  run;

/* Calculate mean and SD for RPKM and transpose*/
proc means data=prep_means noprint;
  var rpkm;
  by fusion_id;
  class trt;
  output out=rpkm_mean mean=mean stddev=SD;
  run;

data rpkm_mean;
  set rpkm_mean;
  where _type_ = 1;
  run;

proc transpose data=rpkm_mean out=flip_rpkm_mean suffix=_mean_rpkm;
  by fusion_id;
  var mean;
  id trt;
  run;

proc transpose data=rpkm_mean out=flip_rpkm_sd suffix=_SD_rpkm;
    by fusion_id;
    var SD;
    id trt;
    run;

/* merge the datasets*/
proc sort data=flip_logrpkm_mean;
    by fusion_id;
    run;

proc sort data=flip_logrpkm_sd;
    by fusion_id;
    run;

proc sort data=flip_rpkm_mean;
    by fusion_id;
    run;

proc sort data=flip_rpkm_sd;
    by fusion_id;
    run;


data ribo.exp_means;
retain fusion_id 
IPmale_SD_logrpkm
IPfemale_SD_logrpkm
InputMale_SD_logrpkm
InputFemale_SD_logrpkm
IPmale_SD_rpkm
IPfemale_SD_rpkm
InputMale_SD_rpkm
InputFemale_SD_rpkm
IPmale_mean_logrpkm
IPfemale_mean_logrpkm
InputMale_mean_logrpkm
InputFemale_mean_logrpkm
IPmale_mean_rpkm
IPFemale_mean_rpkm
InputMale_mean_rpkm
InputFemale_mean_rpkm ;

merge flip_logrpkm_mean flip_logrpkm_sd flip_rpkm_mean flip_rpkm_sd;
  by fusion_id;
  drop _name_;
  run;


data ribo.logrpkm_means;
retain fusion_id
IPmale_mean_logrpkm
IPfemale_mean_logrpkm
InputMale_mean_logrpkm
InputFemale_mean_logrpkm
IPmale_mean_rpkm
IPFemale_mean_rpkm
InputMale_mean_rpkm
InputFemale_mean_rpkm 
;
merge flip_logrpkm_mean flip_rpkm_mean;
by fusion_id;
drop _name_;
run;
