libname fru "!HOME/mclab/Fru_network/sasdata";

proc export data= fru.spike_stack_avg 
            outfile= "!HOME/mclab/Fru_network/output_data/spike_stack_avg.csv" 
            dbms=csv replace;
     putnames=yes;
run;
