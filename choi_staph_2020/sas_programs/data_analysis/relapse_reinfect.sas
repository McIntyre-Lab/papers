
libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";

/*analyze relapse vs reinfection*/

data r_sab_69_unique;
set choi.choi_756;
where sab_status="R";
run;

proc sort data=choi.choi_756;
by sab_status;

proc univariate data=choi.choi_756 normal plot;
by sab_status;
var LOS;
run;

proc univariate data=r_sab_69_unique;
var days;
run;

data relapse_reinfect;
set r_sab_69_unique;
where sab_status="R";
months=days/30;
/*
if months le 1 then class_month=1;
else if months le 2 then class_month=2;
else if months le 3 then class_month=3;
else if months le 4 then class_month=4;
else  class_month=6;

*/

if months le 2 then class_month=2;
else if months le 4 then class_month=4;
else  class_month=6;
if days ge 150 then days_ge_150=1;
else days_ge_150=0;
run;

/*baseline*/



proc freq data=relapse_reinfect;
tables 
pfge
race
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb                
surg30d            
Persistent              
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept;  
 
run;


proc univariate data=relapse_reinfect;
var days;
run;



proc sort data=relapse_reinfect;
by days;
run;

title "R-Sab episode 1-2 month is 30 days";
title1 "5 is anything greater than 150";
title2 "n=69";

proc sort data=relapse_reinfect;
by class_month;
run;

proc freq data=relapse_reinfect;
tables class_month*(
pfge
race
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb                
surg30d            
Persistent              
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept);  
 
run;

proc univariate data=relapse_reinfect normal plot;
by class_month;
var age ;
run;

proc freq data=relapse_reinfect;
tables class_month*outcome;
run;



proc univariate data=relapse_reinfect normal plot;
by class_month;
var apache_total ;
run;


proc freq data=relapse_reinfect;
tables days_ge_150*(
pfge
race
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb                
surg30d            
Persistent              
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept);  
 
run;

proc sort data=relapse_reinfect;
by days_ge_150;
proc univariate data=relapse_reinfect normal plot;
by days_ge_150;
var age ;
run;

proc freq data=relapse_reinfect;
tables class_month*outcome;
run;

proc univariate data=relapse_reinfect normal plot;
by class_month;
var apache_total ;
run;



proc freq data=relapse_reinfect;
tables pfge*class_month*(

dm                      
dialdep                 
neoplasm                
foreignb                
surg30d            
Persistent              
mrsa                    
metendo );  
 
run;
;
run;
