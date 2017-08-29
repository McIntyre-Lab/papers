/* Calculate average coverage for ref hmz, htz, alt hmz */

libname eqtl '/scratch/lfs/patcon/jnewman/sas_analysis/eqtls';


%macro calc_means(setnum);
proc sort data=eqtl.eqtl_data_for_models_&setnum. out=coverage_per_genotype;
   by gene_id feature_id snp_id cell_type genotype;
run;

data coverage_per_genotype2;
   set coverage_per_genotype;
   if log_measurement ne . then measurement=exp(log_measurement)-1;
   else measurement=.;
run;


proc means data=coverage_per_genotype2 noprint;
   by gene_id feature_id snp_id cell_type genotype;
   var measurement;
   output out=avg_measure_per_geno mean=mean stddev=sd;
run;

data geno1 geno2 geno3;
  set avg_measure_per_geno;
  if genotype=0 then output geno1;
  else if genotype=1 then output geno2;
  else if genotype=2 then output geno3;
  keep gene_id feature_id snp_id cell_type mean sd;
run;

data geno1_2;
   format avg_temp f8.3;
   format sd_temp f8.3;
   set geno1;
   length hmz1_expression $30.;
   avg_temp=strip(put(mean, best8.3));
   sd_temp=strip(put(sd, best8.3));
   hmz1_expression=catt(avg_temp, ' ± ',sd_temp);
   rename mean=hmz1_mean sd=hmz1_sd;
   drop avg_temp sd_temp;
run;

data geno2_2;
   format avg_temp f8.3;
   format sd_temp f8.3;
   set geno2;
   length htz_expression $30.;
   avg_temp=strip(put(mean, best8.3));
   sd_temp=strip(put(sd, best8.3));
   htz_expression=catt(avg_temp, ' ± ',sd_temp);
   rename mean=htz_mean sd=htz_sd;
   drop avg_temp sd_temp;
run;

data geno3_2;
   format avg_temp f8.3;
   format sd_temp f8.3;
   set geno3;
   length hmz2_expression $30.;
   avg_temp=strip(put(mean, best8.3));
   sd_temp=strip(put(sd, best8.3));
   hmz2_expression=catt(avg_temp, ' ± ',sd_temp);
   rename mean=hmz2_mean sd=hmz2_sd;
   drop avg_temp sd_temp;
run;
    
proc sort data=geno1_2;
   by gene_id feature_id snp_id cell_type;
proc sort data=geno2_2;
   by gene_id feature_id snp_id cell_type;
proc sort data=geno3_2;
   by gene_id feature_id snp_id cell_type;
run;

data mean_coverage_per_genotype;
   merge geno1_2 geno2_2 geno3_2;
   by gene_id feature_id snp_id cell_type;
run;

data mean_coverage_per_genotype_&setnum.;
   set mean_coverage_per_genotype;
   htz_hmz1_ratio=htz_mean/hmz1_mean;
   hmz2_hmz1_ratio=hmz2_mean/hmz1_mean;
   *drop hmz1_mean hmz1_sd htz_mean htz_sd hmz2_mean hmz2_sd;
run;
%mend;

%calc_means(1);
%calc_means(2);
%calc_means(3);
%calc_means(4);
%calc_means(5);
%calc_means(6);
%calc_means(7);
%calc_means(8);
%calc_means(9);
%calc_means(10);
%calc_means(11);
%calc_means(12);
%calc_means(13);
%calc_means(14);
%calc_means(15);
%calc_means(16);
%calc_means(17);
%calc_means(18);
%calc_means(19);
%calc_means(20);
%calc_means(21);
%calc_means(22);
%calc_means(23);
%calc_means(24);
%calc_means(25);
%calc_means(26);

data eqtl.mean_coverage_per_gt;
   set mean_coverage_per_genotype_1 mean_coverage_per_genotype_2 mean_coverage_per_genotype_3
   mean_coverage_per_genotype_4 mean_coverage_per_genotype_5 mean_coverage_per_genotype_6
   mean_coverage_per_genotype_7 mean_coverage_per_genotype_8 mean_coverage_per_genotype_9
  mean_coverage_per_genotype_10 mean_coverage_per_genotype_11 mean_coverage_per_genotype_12
  mean_coverage_per_genotype_13 mean_coverage_per_genotype_14 mean_coverage_per_genotype_15
  mean_coverage_per_genotype_16 mean_coverage_per_genotype_17 mean_coverage_per_genotype_18
  mean_coverage_per_genotype_19 mean_coverage_per_genotype_20 mean_coverage_per_genotype_21
  mean_coverage_per_genotype_22 mean_coverage_per_genotype_23 mean_coverage_per_genotype_24
  mean_coverage_per_genotype_25 mean_coverage_per_genotype_26;
run;


