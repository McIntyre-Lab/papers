libname cegs "McLab/cegs_ase_paper/sas_data";

data var;
set cegs.compare_ct_var;
run;

data cis_var;
set var;
type="cis";
rename cis_var = variance;
keep fusion_id mating_status cis_var type;
run;

data trans_var;
set var;
type="trans";
rename trans_var = variance;
keep fusion_id mating_status trans_var type;
run;

data var2;
length type $ 6 ;
set cis_var trans_var;
run;

*for R plots;
proc export data=var
outfile = "McLab/cegs_ase_paper/output/cis_trans_variance_data.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=var2
outfile = "McLab/cegs_ase_paper/output/cis_trans_variance_data_stack.csv"
dbms=csv replace;
putnames=yes;
run;

data r2;
set cegs.r2_ct_models;
run;

data cis_r2;
set r2;
type="R2_cis";
rename R2_cis = R2;
keep fusion_id mating_status R2_cis type;
run;

data diff_int_r2;
set r2;
type="R2_diff_int";
rename R2_diff_int = R2;
keep fusion_id mating_status R2_diff_int type;
run;

data diff_trans_r2;
set r2;
type="R2_diff_trans";
rename R2_diff_trans=R2;
keep fusion_id mating_status R2_diff_trans type;
run;

data r2_2;
length type $ 15 ;
set cis_r2 diff_int_r2 diff_trans_r2;
run;


*for R plots;
proc export data=r2
outfile = "McLab/cegs_ase_paper/output/cis_trans_r2_data.csv"
dbms=csv replace;
putnames=yes;
run;

proc export data=r2_2
outfile = "McLab/cegs_ase_paper/output/cis_trans_r2_data_stack.csv"
dbms=csv replace;
putnames=yes;
run;



proc sort data=var;
by mating_status;
run;

*cis and trans variance distribution;
ods listing style=analysis;
goptions reset=all gsfname=graph gsfmode=replace device=gif;
title 'Variance distribution';
proc sgplot data=var;
by mating_status;
density cis_var / legendlabel='cis variance' lineattrs=(pattern=solid);
density trans_var / legendlabel='trans variance' lineattrs=(pattern=2);
keylegend / location=inside position=topright across=1;
xaxis display=(nolabel);
run; quit;


proc sort data=r2;
by mating_status;
run;

*distribution of R2 values;
ods listing style=analysis;
goptions reset=all gsfname=graph gsfmode=replace device=gif;
title 'R2 of cis, trans, interaction distributions';
proc sgplot data=r2;
by mating_status;
density R2_cis / legendlabel='R2 of cis' lineattrs=(pattern=solid);
density R2_diff_int / legendlabel='R2 of diff interaction' lineattrs=(pattern=2);
density R2_diff_trans / legendlabel='R2 of diff trans' lineattrs=(pattern=3);
keylegend / location=inside position=topright across=1;
xaxis display=(nolabel);
run; quit;

*pull out top values for each R2 to make composite graph of ai vs cis, ai vs trans, ai vs c*t;
proc sort data=r2;
by R2_cis;
run; *fusion S16435_SI virgin 0.9948, 0.0022, 0.0001;

proc sort data=r2;
by R2_diff_int;
run; *fusion F587_SI mated 0.0511, 0.7246, -0.0048;

proc sort data=r2;
by R2_diff_trans;
run; *fusion F50887_SI mated -0.0226, 0.0287, 0.9839;




