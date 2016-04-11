/* Normalize the data to centered values based on median */

libname ribo "!MCLAB/arbeitman/arbeitman_ribotag/sas_data";



proc sort data=ribo.all_coverage_counts_with_key2;
  by fusion_id;
  run;

proc sort data= ribo.on_calls_gt_apn0;
  by fusion_id;
  run;

data ribo_all;
  merge ribo.all_coverage_counts_with_key2 (in=in1) ribo.on_calls_gt_apn0 (in=in2);
  by fusion_id;
  if in1;
  run;

data ribo.fusions_all_on2;
  set ribo_all;
  if flag_fusion_all_on0=1;
  logapn=log(apn+1);  
run;

proc sort data=ribo.fusions_all_on2;
  by sample_id;
  run;

proc means data=ribo.fusions_all_on2 noprint;
  by sample_id;
  output out=means
  median(logrpkm)=median_logrpkm
  median(logapn)=median_logapn  
;
  run;

data alldata;
  merge ribo.fusions_all_on2 (in=in1) means (in=in2);
  by sample_id;
  if in1;
  drop _type_ _freq_;
  run;


/* Calculate centered values */
data ribo.norm_median_center2;
  set alldata;
  median_logrpkm_center= logrpkm - median_logrpkm;
  median_logapn_center=logapn - median_logapn;
run;

proc freq data= ribo.norm_median_center2;
tables trt;
run;
