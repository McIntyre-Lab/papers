/* Make SNP-gene index and prepare genotype data for HPC */

/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';


/* Make SNP-gene index */

data snp_index;
   set eqtl.snp_data_w_info;
   keep gene_id snp_id;
run;

proc sort data=snp_index nodup;
   by gene_id snp_id;
run;

data eqtl.snp2gene_index;
   set snp_index;
run;

/* 12762 SNPxGene combinations */


/* Prepare genotype data */

data genotype_data;
  set eqtl.snp_data_w_info;
  keep gene_id snp_id subject_id genotype;
run;

proc sort data=genotype_data;
   by gene_id snp_id;
run;


data eqtl.genotype_data;
   set genotype_data;
run;

