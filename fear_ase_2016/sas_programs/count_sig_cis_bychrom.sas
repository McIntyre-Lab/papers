
libname dmel "S:\McIntyre_Lab\useful_dmel_data\flybase551\sasdata";

proc sort data=cegs.fdr_flag_cis;
by fusion_id;

proc sort data=dmel.Fb551_si_fusions_unique_flagged;
by fusion_id;

*count direction of AI within a line mating status;
data cegs.tests_anno;
merge cegs.fdr_flag_cis (in=in1) dmel.Fb551_si_fusions_unique_flagged(in=in2);
by fusion_id;
if in1;

run;

proc sort data= cegs.tests_anno;
by  mating_status;

proc freq data=cegs.tests_anno noprint;
where flag_fdr_sig=1;
by mating_status;
tables symbol_cat*effect/out=count_sig;
run;
ods results;

proc freq data=count_sig ;
tables effect;
run;

proc sort data=count_sig;
by symbol_cat;
data cis;
set count_sig;
where effect="cis_i";

data trans;
set count_sig;
where effect="trans_i";

data int;
set count_sig;
where effect="cis_i*trans_i";

data all both cis_int trans_int cis2 trans2 int2 oops;
merge cis (in=in1) trans(in=in2) int(in=in3);
by symbol_cat;
if in1 and in2 and in3 then output all;
else if in1 and in2 then output both;
else if in1 and in3 then output cis_int;
else if in2 and in3 then output trans_int;
else if in1 then output cis2;
else if in2 then output trans2;
else if in3 then output int2;
else output oops;
run;

*count both cis and trans;
proc freq data=count_sig;
where effect ne "cis_i*trans_i";
tables symbol_cat;
run;

