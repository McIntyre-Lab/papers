

/*table 3*/

title1 "checking variables";
proc freq data=choi_756;
tables RIDOM_CC*(flag_ridom_cc2 flag_ridom_cc8 ridom_cc2_vs_8) ;
tables flag_ridom_cc2*flag_ridom_cc8;
tables ridom_cc2_vs_8*(flag_ridom_cc2 flag_ridom_cc8);
run;

title1 "S-SAB vs R-SAB";
proc freq data=choi_756;
tables sab_status*(flag_ridom_cc2 flag_ridom_cc8 ridom_cc2_vs_8)/exact ;
run;


title1 "S-SAB vs R-SAB reinfection";
proc freq data=choi_756;
where sab_status="S" or reinfect=1;
tables sab_status*(flag_ridom_cc2 flag_ridom_cc8 ridom_cc2_vs_8)/exact ;
run;

title1 "S-SAB vs R-SAB relapse";
proc freq data=choi_756;
where sab_status="S" or reinfect=0;
tables sab_status*(flag_ridom_cc2 flag_ridom_cc8 ridom_cc2_vs_8)/exact ;
run;


title1 "S-SAB vs R-SAB relapse";
proc freq data=choi_756;
where ridom_cc2_vs_8 ne 0;
tables sab_status*(ridom_cc2_vs_8)/exact ;
run;
