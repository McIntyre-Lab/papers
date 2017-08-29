libname con '!PATCON/sas_data';

/* Make table of demographics */

*MACS data;

proc import datafile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/T1DGC.2011.03_Phenotypic_Current/Analyses Files/MACS_sorting_data.csv"
out=sample_info dbms=csv replace; guessingrows=300;
run;

* Phenotypic data;
proc import datafile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_splicing/T1DGC.2011.03_Phenotypic_Current/Analyses Files/T1DGC.2011.03_Phenotypic_Current_Forms.csv"
out=pheno_data dbms=csv replace; guessingrows=16500;
run;


* Possible covariates;

proc import datafile='/home/jrbnewman/concannon/design_files/diabetes_all_possible_factors.csv'
    out=possible_covariates
    dbms=csv
    replace;
    getnames=yes;
    guessingrows=300;
run;


data design;
  set con.design_by_subject_new;
if name =  '2009-PC-0221' then delete; *sample 75 cd8;
if name =  '2009-PC-0144' then delete; *sample 48 cd4;
if name =  '2009-PC-0236' then delete; *sample 80;
if name =  '2009-PC-0237' then delete; *sample 80;
if name =  '2009-PC-0235' then delete; *sample 80;
run;


data macs;
  set sample_info;
  rename Samp__ID = subject_id;
  drop box;
run;


data pheno_data2;
  set pheno_data;
  keep famid analytic_id psid proband onset duration agee age race1 diabetes;
run;

data covar;
  set possible_covariates;
  keep Name subject_id cell_Type purity live_count dead_count cell_cnt_1e6 total_cells viability;
run;

proc sort data=covar;
  by name;
proc sort data=design;
  by name;
run;

data design_covar;
  merge design (in=in1) covar;
  by name;
  if in1;
run;

proc sort data=design_covar;
  by subject_id;
proc sort data=macs;
  by subject_id;
run;

data design_covar_macs;
  merge design_covar (in=in1) macs (in=in2);
  by subject_id;
  if in1;
run;

proc sort data=design_covar_macs;
  by analytic_id;
proc sort data=pheno_data2;
  by analytic_id;
run;

data design_covar_macs_pheno;
  merge design_covar_macs (in=in1) pheno_data2 (in=in2);
  by analytic_id;
  if in1 and in2;
run;

/* Distribution of:
  male/female, Onset, Age, Duration
  by subject_id

  cell counts by subject_id*cell_type
*/

* subset for by-subject distribs ;

data subject_data;
  set design_covar_macs_pheno;
  keep subject_id age agee sex duration onset race1;
run;

proc sort data=subject_data nodup;
  by subject_id;
run;

*Male female ratio;
proc freq data=subject_data noprint;
  tables sex / out=sex_counts;
run;

* Distribution of age, duration, onset;
proc means data=subject_data noprint;
  var age agee duration onset;
  output out=distrib_by_subject
         mean(age)=mean_age
         stddev(age)=sd_age
         mean(duration)=mean_duration
         stddev(duration)=sd_duration
         mean(onset)=mean_onset
         stddev(onset)=sd_onset;
run;

* Cell counts by cell type and subject;
proc sort data=design_covar_macs_pheno;
   by cell_type;
proc means data=design_covar_macs_pheno noprint;
  by cell_type;
  var live_count dead_count cell_cnt_1e6 total_cells viability;
  output out=distrib_by_cell
         mean(live_count)=mean_live
         stddev(live_count)=sd_live
         mean(dead_count)=mean_dead
         stddev(dead_count)=sd_dead
         mean(viability)=mean_viability
         stddev(viability)=sd_viability
         mean(total_cells)=mean_total
         stddev(total_cells)=sd_total
         mean(cell_cnt_1e6)=mean_count
         stddev(cell_cnt_1e6)=sd_count;
run;

/* Format distribution table */

*sex;
data males females;
  set sex_counts;
  if sex=1 then output males;
  else output females;
  drop percent;
run;

data males2;
  set males;
  rename sex=sex_male count=num_males;
run;

data females2;
  set females;
  rename sex=sex_female count=num_females;
run;

data sex_distrib;
  length factor $255.;
  length distribution $255.;
  merge males2 females2;
  factor="Gender";
  distribution=cat(strip(put(num_males, 10.0)), " males / ",strip(put(num_females,10.0))," females");
  keep factor distribution;
run;

* age, onset, duration;
data age;
  length factor $255.;
  length distribution $255.;
  set distrib_by_subject;
  factor="Age";
  distribution=cat(strip(put(mean_age, 10.0)), " ± " , strip(put(sd_age, 10.0)), " years");
  keep factor distribution;
run;

data onset;
  length factor $255.;
  length distribution $255.;
  set distrib_by_subject;
  factor="Onset age";
  distribution=cat(strip(put(mean_onset, 10.0)), " ± " , strip(put(sd_onset, 10.0)), " years");
  keep factor distribution;
run;

data duration;
  length factor $255.;
  length distribution $255.;
  set distrib_by_subject;
  factor="Disease duration";
  distribution=cat(strip(put(mean_duration, 10.0)), " ± " , strip(put(sd_duration, 10.0)), " years");
  keep factor distribution;
run;

/* Counts by cell */

data total_cells;
  length factor $255.;
  length distribution $255.;
  set distrib_by_cell;
  factor=cat("Total cell count (×10^6) (",strip(cell_type),"+)");
  distribution=cat(strip(put(mean_total, 10.3)), " ± " , strip(put(sd_total, 10.3)));
  if cell_type='' then delete;
  keep factor distribution;
run;

/* Table */

data distrib_table;
  set sex_distrib age onset duration total_cells;
run;

proc export data=distrib_table outfile="!MCLAB/jrbnewman/manuscripts/Newman_T1D_Splicing/reviewer_responses/distribution_table.csv"
     dbms=csv replace;
run;

  
