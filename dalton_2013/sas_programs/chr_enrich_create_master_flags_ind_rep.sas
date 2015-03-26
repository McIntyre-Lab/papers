/* Pull Male Flags */

data male;
    set FRU.Flag_x_induced_repressed_male;
    rename flag_a_ind = flag_male_a_ind;
    rename flag_b_ind = flag_male_b_ind;
    rename flag_c_ind = flag_male_c_ind;
    rename flag_a_rep = flag_male_a_rep;
    rename flag_b_rep = flag_male_b_rep;
    rename flag_c_rep = flag_male_c_rep;
    rename flag_induced = flag_male_ind;
    rename flag_repressed = flag_male_rep;
    keep primary_fbgn flag_a_ind flag_b_ind flag_c_ind flag_a_rep flag_b_rep
         flag_c_rep flag_induced flag_repressed;
    run;

data female;
    set FRU.Flag_x_induced_repressed_female;
    rename flag_a_ind = flag_female_a_ind;
    rename flag_b_ind = flag_female_b_ind;
    rename flag_c_ind = flag_female_c_ind;
    rename flag_a_rep = flag_female_a_rep;
    rename flag_b_rep = flag_female_b_rep;
    rename flag_c_rep = flag_female_c_rep;
    rename flag_induced = flag_female_ind;
    rename flag_repressed = flag_female_rep;
    keep primary_fbgn flag_a_ind flag_b_ind flag_c_ind flag_a_rep flag_b_rep
         flag_c_rep flag_induced flag_repressed;
    run;

data null;
    set FRU.Flag_x_ind_rep_null_male;
    rename flag_induced = flag_null_ind;
    rename flag_repressed = flag_null_rep;
    keep primary_fbgn flag_induced flag_repressed;
    run;

proc sort data=male;
    by primary_fbgn;
    run;

proc sort data=female;
    by primary_fbgn;
    run;

proc sort data=null;
    by primary_fbgn;
    run;


data FRU.flag_ind_rep;
    merge male female null;
    by primary_fbgn;

    if flag_male_a_ind = 1 then flag_male_a = 1;
    else if flag_male_a_rep = 1 then flag_male_a = -1;
    else flag_male_a = 0;

    if flag_male_b_ind = 1 then flag_male_b = 1;
    else if flag_male_b_rep = 1 then flag_male_b = -1;
    else flag_male_b = 0;

    if flag_male_c_ind = 1 then flag_male_c = 1;
    else if flag_male_c_rep = 1 then flag_male_c = -1;
    else flag_male_c = 0;


    if flag_female_a_ind = 1 then flag_female_a = 1;
    else if flag_female_a_rep = 1 then flag_female_a = -1;
    else flag_female_a = 0;

    if flag_female_b_ind = 1 then flag_female_b = 1;
    else if flag_female_b_rep = 1 then flag_female_b = -1;
    else flag_female_b = 0;

    if flag_female_c_ind = 1 then flag_female_c = 1;
    else if flag_female_c_rep = 1 then flag_female_c = -1;
    else flag_female_c = 0;

    if flag_null_ind = 1 then flag_null = 1;
    else if flag_null_rep = 1 then flag_null = -1;
    else flag_null = 0;
    
    run;
