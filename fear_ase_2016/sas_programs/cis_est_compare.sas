
data test4;
set test2;
rename cis_tester=total_cis_tester
cis_i=total_cis_line
trans_t=total_trans_tester
trans_i=total_trans_line;
run;


data test1;
set test4;
Exp=(line+tester)/both;
diff=(tester-line)/both;
tester1=tester/both;
line1=line/both;
*keep geno ms fusionid both line tester exp diff tester1 line1;
run;

proc means data=test1 mean;
by fusionid ms;
var line1 tester1 diff exp;
output out=mean_t1;
run;
data mean_counts;
set mean_t1;
where _stat_="MEAN";
rename _freq_=numobs line1=mean_line1 tester1=mean_tester1 diff=cis_tester exp=mean_exp;
drop _type_;
run;

data mean_counts2;
set mean_counts;
trans_t=2*((mean_tester1)-(mean_exp/2)-cis_tester);
run;

data test2;
merge test1 mean_counts2;
by fusionid ms;

cis_i=((line-tester)/both)+cis_tester;
trans_i=2*(line1-cis_i -(mean_exp/2))-trans_t;
bias=line1/exp;
check_i=(mean_exp/2)+cis_i+((trans_i+trans_t)/2);

run;
*(tester1+line1)/2 is 4.69507755 average for 30 is pop mean;

proc gplot data=test2;

plot bias*cis_i;
plot cis_i*trans_i;
plot line1*cis_i;
plot line1*trans_i;
plot line1*check_i;
plot old_cis_line*cis_i;
plot old_trans_line*trans_i;
plot cis_i*total_cis_line;
plot trans_t*total_trans_tester;
plot trans_i*total_trans_line;

run;

proc means data=test2;
run;


proc univariate data=test2 normal plot;
var cis_i trans_i;
run;
