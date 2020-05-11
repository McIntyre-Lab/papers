
libname choi "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata\data_analysis";

data choi_756;
set choi.choi_756;
run;


data rsab_s_vs_m;
set choi_756;
where sab_status="R";
if num_pairs>1 then flag_num_pairs_ge1=1;
else flag_num_pairs_ge1=0;
run;

proc means data=rsab_s_vs_m p25 median p75;
class flag_num_pairs_ge1;
var age apache_total;
run;


proc npar1way data=rsab_s_vs_m  ;
class flag_num_pairs_ge1;
var age;
run;

proc npar1way data=rsab_s_vs_m  ;
class flag_num_pairs_ge1;
var apache_total;
run;

proc freq data=rsab_s_vs_m ;
tables flag_num_pairs_ge1*(
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

