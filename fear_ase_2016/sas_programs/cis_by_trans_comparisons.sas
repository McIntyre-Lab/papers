filename mymacros "Z:/Mcintyre_Lab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "Z:/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "Z:/useful_dmel_data/flybase551/sasdata";
libname cegs 'Z:/cegs_ase_paper/sas_data';

filename mymacros "McLab/maize_ozone/2014/sas_analysis/macros";
options SASAUTOS=(sasautos mymacros);
%include "McLab/maize_ozone/2014/sas_analysis/macros/iterdataset.sas";
libname dmel "McLab/useful_dmel_data/flybase551/sasdata";
libname cegs 'McLab/cegs_ase_paper/sas_data';

proc sort data=cegs.cis_calls;
by fusion_id;
proc sort data=dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
run;

data cis_calls_w_gene;
merge cegs.cis_calls (in=in1) dmel.fb551_si_fusion_2_gene_id;
by fusion_id;
if in1;
run;

data cis_calls_w_gene_clean;
set cis_calls_w_gene;
if chrom = "2RHet" or chrom= "4" or chrom = "dmel_mitochondrion_genome" then delete;
run;

*plots of trans vs cis by fusion were made using only significant ;
*need to filter: minimum of 5 lines-- run regression-- find slope-- flag is diff than 0 slope-- flag if +/- cis or trans;
data cis_calls_m;
set cis_calls_w_gene_clean;
if mating_status="M";
run;

data cis_calls_V;
set cis_calls_w_gene_clean;
if mating_status="V";
run;

proc freq data=cis_calls_M;
tables fusion_id / out=count_fusions_M;
run;

data less_10;
set count_fusions_M;
if count < 10;
run; *no fusion has < 10 lines significant;

*it looks like there is only a few data points when we plot this so we're going to check it out;
data chk_F20713;
set cis_calls_M;
if fusion_id = "F20713_SI";
run;

proc gplot data=chk_F20713;
plot trans_line*cis_line;
run; quit; *the R plots are not quite correct, we don't need to drop anything;

*check V;
proc freq data=cis_calls_V;
tables fusion_id / out=count_fusions_V;
run;

data less_10;
set count_fusions_V;
if count < 10;
run; *none < 10;

/* no filter for mated or virgin */

proc sort data=cis_calls_M;
by fusion_id;
proc sort data=cis_calls_V;
by fusion_id;
run;

data test_M;
set cis_calls_M;
if fusion_id = "F10001_SI" or fusion_id = "S47599_SI";
run;

proc reg data=test_M OUTEST=reg_output ALL;
by fusion_id;
      model cis_line = trans_line;
	  ods output ANOVA=ANOVA CORR=CORR FITSTATISTICS=FITSTAT;
*=fitstatistics out=ANOVA;
   run; 

proc reg data=cis_calls_M ;
by fusion_id;
      model cis_line = trans_line;
	  ods output ANOVA=ANOVA_M CORR=CORR_M ;
   run; 

proc reg data=cis_calls_V ;
by fusion_id;
      model cis_line = trans_line;
	  ods output ANOVA=ANOVA_V CORR=CORR_V ;
   run; 


