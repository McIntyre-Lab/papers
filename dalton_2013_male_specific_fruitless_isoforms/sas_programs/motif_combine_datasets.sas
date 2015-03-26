libname '!MCLAB/arbeitman_fru_network/sasdata';

data female;
    set fru.flag_x_induced_repressed_female;
    rename flag_a_ind     = female_a_ind;
    rename flag_b_ind     = female_b_ind;
    rename flag_c_ind     = female_c_ind;
    rename flag_a_rep     = female_a_rep;
    rename flag_b_rep     = female_b_rep;
    rename flag_c_rep     = female_c_rep;
    rename flag_induced   = female_ind;
    rename flag_repressed = female_rep;
    keep primary_fbgn symbol flag_a_ind flag_b_ind flag_c_ind flag_a_rep flag_b_rep flag_c_rep flag_induced flag_repressed;
    run;

data male;
    set fru.flag_x_induced_repressed_male;
    rename flag_a_ind     = male_a_ind;
    rename flag_b_ind     = male_b_ind;
    rename flag_c_ind     = male_c_ind;
    rename flag_a_rep     = male_a_rep;
    rename flag_b_rep     = male_b_rep;
    rename flag_c_rep     = male_c_rep;
    rename flag_induced   = male_ind;
    rename flag_repressed = male_rep;
    keep primary_fbgn symbol flag_a_ind flag_b_ind flag_c_ind flag_a_rep flag_b_rep flag_c_rep flag_induced flag_repressed;
    run;

data null_male;
    set fru.flag_x_ind_rep_null_male;
    rename flag_induced   = null_ind;
    rename flag_repressed = null_rep;
    keep primary_fbgn symbol flag_induced flag_repressed;
    run;

proc sort data=female;
    by primary_fbgn;
    run;
proc sort data=male;
    by primary_fbgn;
    run;
proc sort data=null_male;
    by primary_fbgn;
    run;

data fru.all_motif_flags;
    merge female male null_male;
    by primary_fbgn;
    run;

proc export data=fru.all_motif_flags
            outfile='!MCLAB/arbeitman_fru_network/reports_external/all_motif_flags.csv'
            dbms=csv replace;
            putnames = yes;
            run;

/* example subset */

data mysub;
    set fru.all_motif_flags;
    if male_a_ind = 1 and male_b_ind = 1 and null_ind = 1 and male_c_ind = 0;
    run;
