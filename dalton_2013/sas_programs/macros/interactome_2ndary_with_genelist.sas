

*macro for merging a genelist to fly secondary interactome ;
*need prep genelist:
	genelist_int:
		rename symbol = int_gene and rename significance flag = flag_sig_int 
		keep flag_sig and symbol
	genelist_focal:
		rename symbol = focal_gene  
		keep flag_sig and symbol;
		

%macro fly_2ndary_genelist (genelist_int, genelist_focal, suffix, cutoff1, cutoff2, cutoff3) ;

data secondary ;
set int.fly_2ndary_interactome ;
pres_2ndary = 1;
run;

proc sort data = secondary ;
by int_gene ;
proc sort data = &genelist_int ;
by int_gene ;
run;

proc sql ;
  create table interact_sql as
      select *
      from &genelist_int natural full join
          secondary;
  quit;	

data flagged1 not_flagged  ;
length focal_gene $12. ;
set interact_sql;
if flag_sig_int  = 1 or flag_sig_int  = 0 then output flagged1 ;
rename flag_sig_int  = flag_sig_int_&suffix ;
else output not_flagged;
run;

data flagged other ;
set flagged1;
if pres_2ndary = 1  then output flagged ;
else output other ;
run;	

proc sort data = flagged   ;
by focal_gene;
run;

proc means data = flagged noprint ;
by focal_gene ;
var flag_sig_int_&suffix ;
output out= count_net n=num_int_&suffix sum=num_sig_int_&suffix ;
run;


proc sort data = count_net ;
by focal_gene ;
proc sort data = &genelist_focal ;
by focal_gene ;
run ;

data interaction_1 other  ;
merge count_net (in=in1) &genelist_focal (in=in2)  ;
by focal_gene ;
if in1 and in2 then output interaction_1 ;
else output other ;
run ;

data interaction_2 ;
set interaction_1 ;
rename flag_sig = sig_focal_&suffix ;
percent_sig_&suffix =num_sig_int_&suffix/num_int_&suffix*100;
drop _freq_ _type_ ;
run ;

data interaction_2ndary_&suffix ;
set interaction_2;
if  percent_sig_&suffix ge &cutoff1 then &suffix._&cutoff1=1;
	else &suffix._&cutoff1=0;
if percent_sig_&suffix=. then &suffix._&cutoff1=.;

if  percent_sig_&suffix  ge &cutoff2 then &suffix._&cutoff2=1;
	else &suffix._&cutoff2=0;
if percent_sig_&suffix=. then &suffix._&cutoff2=.;

if  percent_sig_&suffix  ge &cutoff3 then &suffix._&cutoff3=1;
	else &suffix._&cutoff3=0;
if percent_sig_&suffix=. then &suffix._&cutoff3=.;

if num_int_&suffix >1 and &suffix._&cutoff1=1 then network_&suffix = 1;
	else network_&suffix=0;
if percent_sig_&suffix=. then network_&suffix=.;
run;

proc sort data = interaction_2ndary_&suffix ;
by focal_gene;
run;

%mend ;



