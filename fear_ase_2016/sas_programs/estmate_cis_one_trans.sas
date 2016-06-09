*libname cegs "S:\SHARE\\McIntyre_Lab\cegs_ase_paper\sas_data";
libname cegs "!MCLAB/cegs_ase_paper/sas_data";

libname dmel "/home/mcintyre/S/SHARE/McIntyre_Lab/useful_dmel_data/flybase551/sasdata";
libname cegs "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/sas_data";


libname gwas '/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_gwas_paper/pipeline_output';


proc contents data=cegs.Clean_ase_stack;
run;

/* this should take into account the qsim*/
data sig;
set cegs.Clean_ase_stack;
if flag_all_ai=1;
run;

proc freq data=sig noprint;
tables fusion_id/out=count_lines;
run;

data keep_for_cis;
set count_lines;
if count ge 10;
run;

proc sort data=keep_for_cis;
by fusion_id;
proc sort data=sig;
by fusion_id;

data cis_data;
*merging all data instead of just significant for ai;
merge cegs.Clean_ase_stack(in=in1) keep_for_cis(in=in2);
by fusion_id;
if in1 and in2;
run;

data cis_data_normalized;
set cis_data;
sum_ase=sum_line+sum_tester;
line1=sum_line/sum_both;
tester1=sum_tester/sum_both;
diff=tester1-line1;
test=line1+tester1;
overall_expression=sum_ase+sum_both;
run;

proc sort data=cis_data_normalized;
by mating_status fusion_id;

proc means data=cis_data_normalized noprint;
by mating_status fusion_id;
var overall_expression;
output out=sum_OE sum=total_express;

proc print data=cis_est_one_trans (obs=10);
run;

proc means data=cis_data_normalized noprint;
class fusion_id mating_status;
var line1 tester1 diff;
output out=sum sum(line1)=pop_sum_line sum(tester1)=pop_sum_tester mean(diff)=C_t;
run;

data pop_mean;
set sum;
where _type_=3;
sum_all=pop_sum_line+(pop_sum_tester/_freq_);
*allele_mean=sum_all/_freq_; *this should be .5 I think;
allele_mean=pop_sum_line/_freq_;
drop _type_;
run;

proc sort data=pop_mean ;
by fusion_id mating_status;
proc sort data=cis_data_normalized ;
by fusion_id mating_status;
proc sort data=sum_oe;
by fusion_id mating_status;

data cis_est_one_trans;
merge cis_data_normalized pop_mean sum_oe;
by fusion_id mating_status;
c_i=(line1-tester1)+C_t;

T_t_1=2*((pop_sum_tester/_freq_)-allele_mean-C_t);
T_i_1a=2*(line1-allele_mean-c_i)-T_t_1;

T_t_2=2*((pop_sum_tester/_freq_)-C_t);
T_i_2a=2*(line1-c_i)-T_t_2;

T_t_3=2*((pop_sum_tester/_freq_)-0.5-C_t);
T_i_3a=2*(line1-0.5-c_i)-T_t_3;

if c_i ge 0 then direction_cis="+";
        else direction_cis="-";
if t_i_1a ge 0 then direction_trans="+";
        else direction_trans="-";

t_pop_i=sum_total-pop_sum_tester-allele_mean-C_t;

run;

*checking;
proc means data =  cis_est_one_trans noprint;
by fusion_id mating_status;
var  line1 tester1 sum_total pop_sum_tester allele_mean c_t t_pop_i;
output out=t_pop mean( line1 tester1 sum_total pop_sum_tester allele_mean c_t t_pop_i ) =/autoname;
run;

proc print data=cis_est_one_trans (obs=10);
run;


proc print data=t_pop  (obs=10);
run;

data check;


proc univariate data=t_pop normal plot;
var t_pop_mean;
run;

proc sort data=t_pop;
by fusion_id mating_status;

data cis_est_trans_pop;
merge cist_est_one_trans t_pop;
by fusion_id mating_status;
run;



proc freq data=cegs.cis_est_v13;
tables direction_cis*direction_trans/agree;
run;

proc gplot data=cegs.cis_est_v13;
plot q5_mean_theta*c_i;
plot line1*c_i;
plot c_i*(t_i_1a t_i_2a t_i_3a);
plot t_i_1a*(t_i_2a t_i_3a);
run;

