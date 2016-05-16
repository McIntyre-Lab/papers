PROC IMPORT OUT= WORK.diverge 
            DATAFILE= "/home/agerken/McLab/alison_g/cegs_ase_explore/RAL
_pi_USC_browser.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;
