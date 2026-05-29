libname seq "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/t1d_case_control_cellType/sasdata" ;

/*
flag short frags -- less than 10 bp
*/


/* import frag annotation flag shrt frags   == 558,791 obs/frags */

proc import datafile = "//TB14/TB14/t1d_case_control_cellType/allChr_isoseq_FSM_ISM_NIC_event_analysis_ef.csv"
out = anno
dbms = csv replace ;
guessingrows = MAX ;
run;


proc contents data = anno ; run;

/* flag short fragments */
data anno2 ;
set anno ;
frag_length = ef_end - ef_start ;
if frag_length < 10 then flag_frag_less_10bp = 1 ;
else flag_frag_less_10bp = 0; 
rename ef_id = featureID ;
keep ef_id flag_frag_less_10bp  ;
run ; 

proc freq data = anno2 ;
tables flag_frag_less_10bp ;
run; /* 25,050 frags out of 151,719 that are less than 10 bp (16.5%) */

data seq.flag_shrt_frags ;
set anno2 ;
run;

proc freq data = seq.flag_shrt_frags ; 
tables featureID / out = cnt_shrts ;
run;


