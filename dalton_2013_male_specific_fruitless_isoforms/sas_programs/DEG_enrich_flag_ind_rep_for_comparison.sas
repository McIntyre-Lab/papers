/******************************************************************************** 
* In order to do the comparison that we are wanting to do, I need to create new
* flags (-1,0,1) where (-1) is repressed, (0) is no change, (1) is induced.
********************************************************************************/

data male;
    set FRU.flag_x_induced_repressed_male;
    if flag_a_ind = 1 then flag_a = 1; 
    else if flag_a_rep = 1 then flag_a = -1;
    else flag_a = 0;
    if flag_b_ind = 1 then flag_b = 1; 
    else if flag_b_rep = 1 then flag_b = -1;
    else flag_b = 0;
    if flag_c_ind = 1 then flag_c = 1; 
    else if flag_c_rep = 1 then flag_c = -1;
    else flag_c = 0;
    keep primary_fbgn symbol flag_a flag_b flag_c;
    run;


data null;
    set FRU.flag_x_ind_rep_null_male;
    if flag_induced = 1 then flag_null = 1; 
    else if flag_repressed = 1 then flag_null = -1;
    else flag_null = 0;
    keep primary_fbgn symbol flag_null;
    run;

proc sort data=male;
    by primary_fbgn;
    run;

proc sort data=null;
    by primary_fbgn;
    run;

data FRU.flag_male_status;
    merge male (in=in1) null(in=in2);
    by primary_fbgn;
    run;
    
proc datasets nolist;
    delete male; 
    delete null;
    run; quit;
