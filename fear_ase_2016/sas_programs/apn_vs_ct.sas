
proc gplot data=cegs.Cis_est_v13;
plot mean_apn*(c_i t_i_1a);
run;


proc gplot data=cegs.Cis_est_v13;
where flag_all_ai=1;
plot mean_apn*(c_i t_i_1a);
run;

proc gplot data=cegs.Cis_est_v13;
where mean_apn le 100;
plot mean_apn*(c_i t_i_1a);
run;


proc gplot data=cegs.Cis_est_v13;
where flag_all_ai=1 and mean_apn le 100;
plot mean_apn*(c_i t_i_1a);
run;
