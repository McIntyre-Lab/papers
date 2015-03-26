libname fru '!MCLAB/Fru_network/sasdata';

proc export data= fru.results_plus_gov2 
            outfile= "/home/jfear/mclab/Fru_network/reports/results_plus_go_v2_20120721.csv" 
            dbms=csv replace;
     putnames=yes;
run;
