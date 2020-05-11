data check_none;
set choi_756;
where sab_status="R";
sum=
dm                      
+dialdep                 
+idu                     
+neoplasm                
+transpat                
+steroid                 
+hiv                     
+foreignb   
+surg30d            
+Persistent   
+mrsa                    
+metendo                 
+metabsc                 
+metarth                 
+metsept;
run;

proc means data=check_none min median max;
var sum;
run; 

data check_none_iso;
set check_none;
where sum=0;
run;


proc export data=check_none_iso
            OUTFILE= "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\output\check_none.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

 

