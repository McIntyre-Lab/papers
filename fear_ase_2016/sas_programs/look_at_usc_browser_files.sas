libname cegs "McLab/alison_g/cegs_ase_explore/sas_data";

PROC IMPORT OUT= WORK.diverge 
            DATAFILE= "/home/agerken/McLab/alison_g/cegs_ase_explore/RAL
_pi_USC_browser.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

PROC IMPORT OUT= WORK.tajima 
            DATAFILE= "/home/agerken/McLab/alison_g/cegs_ase_explore/RAL
_TsD_USC_browser.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data tajima2;
set tajima;
rename name=tajima_d;
run;

data diverge2;
set diverge;
rename name=diverge_pi;
run;
