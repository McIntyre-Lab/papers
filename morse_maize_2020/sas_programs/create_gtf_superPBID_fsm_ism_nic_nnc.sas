libname pb "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata";

/*
subset sqanti post filter gtf file to only fsm, ism, nic and nnc

*/


proc import datafile =  "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/sqanti_filtered_corrected.gtf"
out = gtf 
dbms = tab replace ;
guessingrows = MAX ;
getnames = no ;
run ;


data gtf2 ;
retain pbid ;
set gtf ;
length isoform $13. ;
tr = scan(var9, 2, " ");
isoform = dequote(tr) ;
drop tr ;
run;

data list ;
set pb.fsm_ism_isoform pb.nic_nnc_isoform ;
run;

proc sort data = gtf2 ;
by isoform ;
proc sort data = list ;
by isoform ;
run;

data gtf3 not_on_list;
merge list (in=in1) gtf2 (in=in2) ;
if in1 then output gtf3 ;
else output not_on_list ;
by isoform ;
run;

data subset_gtf ;
set gtf3 ;
drop isoform ;
run ;
    /*  33,267 obs in list
        326,870 obs in gtf file
        subset gtf has 262,280 obs   */

data pb.pbID_fsm_ism_nic_nnc_gtf ;
set subset_gtf ;
run ;

data new_pbID_fsm_ism_nic_nnc_gtf ;
retain var1 var2 var3 var4 var5 newvar6 var7 newvar8 var9 ;
set pb.pbID_fsm_ism_nic_nnc_gtf  ;
if var6 = . then newVar6 = ".";
if var8 = . then newVar8 = ".";
drop var6 var8 ;
run ;

data pb.pbID_fsm_ism_nic_nnc_gtf ;
set new_pbID_fsm_ism_nic_nnc_gtf  ;
run ;

data _null_ ;
file "!MCLAB/maize_ozone_FINAL/2018/PacBio/sqanti_classification_category_subset/pbID_fsm_ism_nic_nnc.gtf" delimiter = '09'x;
set   pb.pbID_fsm_ism_nic_nnc_gtf ;
put (_all_) (+0);
run;



