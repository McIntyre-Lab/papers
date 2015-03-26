libname fru '!HOME/mclab/Fru_network/sasdata';
/* 
* I was provided with ng per 100ng of RNA. For the ERCC spikes we actually use
* a concentration adjusted for length. For now I am just going to use ng
* amount, initially this should be good enough
*/

data fru.spike_design_file;
    input spike_id $ conc ;
    datalines;
    EC3 0.04
    EC15 0.02
    EC17 0.01
    EC13 0.005
    EC18 0.0025

    ;
    run;
