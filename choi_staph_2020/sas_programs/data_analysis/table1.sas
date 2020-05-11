libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";

data choi_756;
set choi.choi_756;
run;

/*Table 1 a*/

title1 "S-SAB vs R-SAB";

proc means data=choi_756  p25 median p75;
class sab_status;
var age apache_total;
run;


proc npar1way data=choi_756  ;
class sab_status;
var age;
run;

proc npar1way data=choi_756  ;
class sab_status;
var apache_total;
run;

proc freq data=choi_756;
tables sab_status*(
female
race_3level
race_2level
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb   
route 
surg30d            
Persistent   
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept)/exact;  
 
run;



title1 "S-SAB vs reinfection";

data ssab_vs_reinfect;
set choi_756;
where sab_status="S" or reinfect=1;
run;


proc means data=ssab_vs_reinfect p25 median p75;
class sab_status;
var age apache_total;
run;


proc npar1way data=ssab_vs_reinfect  ;
class sab_status;
var age;
run;

proc npar1way data=ssab_vs_reinfect ;
class sab_status;
var apache_total;
run;

proc freq data=ssab_vs_reinfect;
tables sab_status*(
female
race_3level
race_2level
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb                
route
surg30d            
Persistent              
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept)/exact;  
 
run;


title1 "S-SAB vs relapse";

data ssab_vs_relapse;
set choi_756;
where sab_status="S" or reinfect=0;
run;


proc means data=ssab_vs_relapse p25 median p75;
class sab_status;
var age apache_total;
run;

proc npar1way data=ssab_vs_relapse  ;
class sab_status;
var age;
run;

proc npar1way data=ssab_vs_relapse;
class sab_status;
var apache_total;
run;

proc freq data=ssab_vs_relapse;
tables sab_status*(
female
race_3level
race_2level
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb
route 
surg30d            
Persistent              
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept)/exact;  
 
run;



title1 "Relapse vs reinfection";

data rsab;
set choi_756;
where sab_status="R";

proc means data=rsab  p25 median p75;
class reinfect;
var age apache_total;
run;

proc npar1way data=rsab  ;
class reinfect;
var age;
run;

proc npar1way data=rsab ;
class reinfect;
var apache_total;
run;


proc freq data=rsab;
tables reinfect*(
female
race_2level
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb
route 
surg30d            
Persistent              
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept
)/exact;  
 
run;

title1 "check if early relapse is important";
title2 " 1 is early relapse less than 60 days 0 is relapse";



proc means data=rsab  p25 median p75;
class early_relapse;
var age apache_total;
run;

proc npar1way data=rsab  ;
class early_relapse;
var age;
run;

proc npar1way data=rsab ;
class early_relapse;
var apache_total;
run;

proc freq data=rsab;
tables early_relapse*(
female
race_2level
dm                      
dialdep                 
idu                     
neoplasm                
transpat                
steroid                 
hiv                     
foreignb
route 
surg30d            
Persistent              
mrsa                    
metendo                 
metabsc                 
metarth                 
metsept
)/exact;  
 
run;
