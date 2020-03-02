
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";

/* up - down patterns --TRANSCRIPTS 
    is pattern consistent across genotypes?  
*/

data trans_cnts3 ;
set tappas.tappas_results_transcripts ;
where regulated_B73  ne "" or regulated_Mo17  ne "" or regulated_C123  ne "" or regulated_Hp301  ne "" or regulated_NC338  ne "" ; 
count = 5 - cmiss(regulated_b73,  regulated_c123, regulated_hp301, regulated_mo17, regulated_nc338) ; 
keep name_description regulated_: count;
run;   /* 4922 genes resulted in 1 or more genotype */

proc sort data = trans_cnts3 ;
by name_description ;
run;

data trans_cnts4 ;
retain pattern ;
sets trans_cnts3 ;
pattern = compress(regulated_b73||'|'||regulated_c123||'|'||regulated_hp301||'|'||regulated_mo17||'|'||regulated_nc338);
run ;

data consistent ;
retain consistent ;
length consistent $3.;
set trans_cnts4 ;
if find(pattern,"DOWN") ge 1 and find(pattern,"UP") ge 1 then consistent = 'no' ;
else consistent = 'yes' ;
run;

proc freq data = consistent ;
tables consistent ;
run;
/*
consistent    Frequency
-----------------------
no                  51
yes               7405
*/



/*
check consistency of detection across genotypes for transcripts 
*/

data trans_cnts3 ;
set tappas.tappas_results_transcripts ;
where dea_result_B73  ne "" or dea_result_Mo17  ne "" or dea_result_C123  ne "" or dea_result_Hp301  ne "" or dea_result_NC338  ne "" ; 
count = 5 - cmiss(dea_result_b73,  dea_result_c123, dea_result_hp301, dea_result_mo17, dea_result_nc338) ; 
keep name_description dea_result_: count;
run;   /* 23,194 detected*/

proc freq data = trans_cnts3 ;
tables count ;
run;/*
count    Frequency     Percent     Frequency
--------------------------------------------
    1        2927       12.62          2927
    2        2299        9.91          5226
    3        2298        9.91          7524
    4        2644       11.40         10168
    5       13026       56.16         23194
*/



proc sort data = trans_cnts3 ;
by name_description ;
run;

data trans_cnts4 ;
retain pattern ;
set trans_cnts3 ;
length b73_on $3. ;
length c123_on $3. ;
length hp301_on $3. ;
length mo17_on $3. ;
length nc338_on $3. ;

if dea_result_b73 ne "" then b73_on = 'on'; else b73_on = 'off';
if dea_result_c123 ne "" then c123_on = 'on'; else c123_on = 'off';
if dea_result_hp301 ne "" then hp301_on = 'on'; else hp301_on = 'off';
if dea_result_mo17 ne "" then mo17_on = 'on'; else mo17_on = 'off';
if dea_result_nc338 ne "" then nc338_on = 'on'; else nc338_on = 'off';
pattern = compress(b73_on ||'|'||c123_on||'|'||hp301_on||'|'||mo17_on||'|'||nc338_on) ;

keep name_description  b73_on c123_on hp301_on mo17_on nc338_on pattern ;
run;


data trans_cnt5 ;
retain pattern  ;
set trans_cnts3 ;
if dea_result_b73 = "Not DE" then b73_trt = 'amb'; else if dea_result_b73 = "DE" then b73_trt = 'oz'; else b73_trt ='off';
if dea_result_c123 = "Not DE" then c123_trt = 'amb'; else if dea_result_c123 = "DE" then c123_trt = 'oz'; else c123_trt ='off';
if dea_result_hp301 = "Not DE" then hp301_trt = 'amb'; else if dea_result_hp301 = "DE" then hp301_trt = 'oz'; else hp301_trt ='off';
if dea_result_mo17 = "Not DE" then mo17_trt = 'amb'; else if dea_result_mo17 = "DE" then mo17_trt = 'oz'; else mo17_trt ='off';
if dea_result_nc338 = "Not DE" then nc338_trt = 'amb'; else if dea_result_nc338 = "DE" then nc338_trt = 'oz'; else nc338_trt ='off';
pattern = compress(b73_trt||'|'||c123_trt||'|'||hp301_trt||'|'||mo17_trt||'|'||nc338_trt) ;
keep name_description b73_trt c123_trt hp301_trt mo17_trt nc338_trt pattern ;
run ;


data consistent ;
retain consistent ;
length consistent $3.;
set trans_cnt5 ;
if find(pattern,"amb") ge 1 and find(pattern,"oz") ge 1 then consistent = 'no' ;
else consistent = 'yes' ;
run;

proc freq data = consistent ;
tables pattern consistent ;
run;
/*
 *** ignoring off, if all amb or all oz then consistent = yes, else if  mix of amb and oz then consistent = no
consistent    Frequency     Percent     Frequency
-------------------------------------------------
no                6616       28.52          6616
yes              16578       71.48         23194
*/



