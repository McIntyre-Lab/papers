
libname cegs "S:\McIntyre_Lab\cegs_ase_paper\sas_data";
libname ase "Y:\ase_extended\sas_data";

proc contents data=cegs.Clean_ase_stack;
run;

proc contents data=ase.gold_key;
run;


proc sort data=cegs.Clean_ase_stack;
by fusion_id line;
run;

proc sort data=ase.gold_key;
by fusion_id line;
run;


data check_bias;
merge cegs.Clean_ase_stack (in=in1) ase.gold_key(in=in2) ;
by fusion_id line;
if in2;
run;

proc sgplot data= check_bias;
where mating_status="M";
density q5_mean_theta;
xaxis grid values =(0 to 1 by .05);
refline .5 / axis=x;
run;
