/* Stack fusion and splicing events into one expression dataset */

libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';


data splicing_counts;
   set eqtl.splicing_counts;
   log_measurement=log(depth+1);
   keep Name cell_type subject_id event_id log_measurement gene_id;
   rename event_id=feature_id;
run;

data fusion_counts;
   set eqtl.fusion_counts_for_eqtls;
   keep Name cell_type subject_id fusion_id log_q3_q3_apn gene_id;
   rename fusion_id=feature_id log_q3_q3_apn=log_measurement;
run;

data splicing_counts2;
   retain Name cell_type subject_id gene_id feature_id log_measurement;
   set splicing_counts;
run;

data fusion_counts2;
   retain Name cell_type subject_id gene_id feature_id log_measurement;
   set fusion_counts;
run;

data expression_data;
    set splicing_counts2 fusion_counts2;
run;

data gene_list;
    set eqtl.snp2gene_index;
    keep gene_id;
run;

proc sort data=gene_list nodup;
   by gene_id;
proc sort data=expression_data;
   by gene_id;
run;

data eqtl.expression_data;
   merge gene_list (in=in1) expression_data;
   by gene_id;
  if in1;
run;

proc sort data=eqtl.expression_data;
   by gene_id subject_id;
run;




