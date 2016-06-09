

data cis_estimates_plus_flags;
set cis_data_estimates;
if q5_mean_theta ge 0.7 or q5_mean_theta le 0.3 then bias_level=3;
else if q5_mean_theta ge 0.6 or q5_mean_theta le 0.4 then bias_level=2;
else bias_level=1;
/* trasn distribution
95% 6.97266E-01 
90% 3.26725E-01 
75% Q3 7.84601E-02 
50% Median -1.50088E-02 
25% Q1 -1.26643E-01 
10% -4.05678E-01 
5% -8.53626E-01 
*/
if trans_i le -4.05678E-01 or trans_i ge 3.26725E-01 then trans_level=3;
else if trans_i le -1.26643E-01 or trans_i ge 7.84601E-02  then trans_level=2;
else trans_level=1;
/* cis distribution
100% Max 1.91613E+02 
99% 2.27601E+00 
95% 3.31084E-01 
90% 1.60126E-01 
75% Q3 4.88047E-02 
50% Median 2.90906E-03 
25% Q1 -3.80524E-02 
10% -1.52432E-01 
5% -3.09231E-01 
1% -1.87041E+00 
0% Min -4.63853E+02 
*/
if cis_i le -1.52432E-01 or cis_i ge 1.60126E-01 then cis_level=3;
else if cis_i le -3.80524E-02 or cis_i ge 4.88047E-02  then cis_level=2;
else cis_level=1;
run;

proc freq data=cis_estimates_plus_flags;
tables cis_level*trans_level;
tables cis_level*bias_level;
tables trans_level*bias_level;
run;

proc sort data=cis_estimates_plus_flags;
by cis_level;
run;

proc univariate data=cis_estimates_plus_flags normal plot;
by cis_level;
var sum_line sum_tester ;
run;


proc sort data=cis_estimates_plus_flags;
by trans_level;

proc univariate data=cis_estimates_plus_flags normal plot;
by trans_level;
var sum_line sum_tester ;
run;


proc sort data=cis_estimates_plus_flags;
by cis_level;
run;

proc univariate data=cis_estimates_plus_flags normal plot;
where sum_line le 1000 and sum_tester le 1000;
by cis_level;
var sum_line sum_tester ;
run;


proc sort data=cis_estimates_plus_flags;
by trans_level;

proc univariate data=cis_estimates_plus_flags normal plot;
where sum_line le 1000 and sum_tester le 1000;
by trans_level;
var sum_line sum_tester ;
run;
