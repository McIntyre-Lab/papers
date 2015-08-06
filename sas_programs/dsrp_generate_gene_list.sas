/********************************************************************************
* Multicollinearity is going to be a large problem. It will be best to
* collapse isoforms if they have a large correlation coefficient (>=0.8).
* 
* This script calculates the correlation between isoforms and flags them if
* they have a correlation >=0.8.
********************************************************************************/

/*
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create Gene list */
proc sort data=SEM.dsrp_sbs_symbol;
    by patRIL matRIL;
    run;

proc transpose data=SEM.dsrp_sbs_symbol out=flip;
    by patRIL matRIL;
    var :;
    run;
    
proc sort data=flip;
    by patRIL matRIL _name_;
    run;

data dsrp_stack;
    set flip;
    keep patRIL matRIL _name_;
    run;

data gene_list;
    format symbol $12.;
    set dsrp_stack;
    count = index(_name_, '_');
    if count >0 then do;
        symbol = trim(substr(_name_, 1, count-1));
    end;
    else symbol = _name_;
    pos=prxmatch("/^fl_2_d/",_name_); 
    if pos > 0 then symbol = 'fl_2_d'; 
    keep symbol ;
    run;

proc sort data=gene_list nodupkey;
    by symbol;
    run; * 7422 obs;

data SEM.dsrp_gene_list;
    set gene_list;
    run;

/* Clean Up */
    proc datasets nolist;
        delete gene_list ;
        delete dsrp_stack;
        delete flip;
        run; quit;
