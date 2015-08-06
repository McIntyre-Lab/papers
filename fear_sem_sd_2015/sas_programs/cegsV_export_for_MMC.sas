/********************************************************************************
* Prepare datset for upload to the MMC website. Genes should be rows and
* observations should be columns.
********************************************************************************/
/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname cegs '!MCLAB/cegs_sergey/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

proc sort data=SEM.cegsV_sex_det_stack;
    by line symbol_cat fusion_id;
    run;

data sex;
    set SEM.cegsV_sex_det_stack;
    if symbol_cat eq 'Spf45' or
    symbol_cat eq 'Sxl' or
    symbol_cat eq 'Yp1' or
    symbol_cat eq 'Yp2' or
    symbol_cat eq 'Yp3' or
    symbol_cat eq 'fl_2_d' or
    symbol_cat eq 'fru' or
    symbol_cat eq 'her' or
    symbol_cat eq 'ix' or
    symbol_cat eq 'snf' or
    symbol_cat eq 'tra' or
    symbol_cat eq 'tra2' or
    symbol_cat eq 'vir' ;
    symbol = trim(symbol_cat) || '_' || trim(fusion_id);
    drop fbgn_cat fusion_id symbol_cat;
    run;

proc sort data=sex;
    by symbol;
    run;

proc transpose data=sex out=flip;
    by symbol;
    var mean_exp;
    id line;
    run;

data flip2;
    set flip;
    drop _name_;
    run;

proc export data=flip2 outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/mmc/cegsV_sex_det_data_for_mmc.csv' dbms=csv replace;
putnames=yes;
run;
