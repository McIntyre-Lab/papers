
*libname fru 'z:\arbeitman\arbeitman_fru_network\sasdata';

data flag_ind_rep;
    set FRU.flag_ind_rep;
    keep primary_fbgn flag_male_a_ind flag_male_b_ind flag_male_c_ind
    flag_male_a_rep flag_male_b_rep flag_male_c_rep flag_female_a_ind
    flag_female_b_ind flag_female_c_ind flag_female_a_rep flag_female_b_rep
    flag_female_c_rep flag_null_ind flag_null_rep ;
    run;

data flag_motif;
    set FRU.motif_flags_and_cnts;
    keep primary_fbgn flag_fru_a_motif flag_fru_b_motif flag_fru_c_motif;
    run;

data region;
    set FRU.motif_search_regions;
    keep primary_fbgn region_length;
    run;

proc sort data=flag_ind_rep;
    by primary_fbgn;
    run;

proc sort data=flag_motif;
    by primary_fbgn;
    run;

proc sort data=region;
    by primary_fbgn;
    run;

data merged;
    merge flag_ind_rep flag_motif region;
    by primary_fbgn;
    run;

data motif;
    set merged;
    if region_length> 10000 then flag_long=1;
    else if region_length<0 then flag_long=.;
    else flag_long=0;
    if flag_fru_a_motif=1 or flag_fru_b_motif=1 or flag_fru_c_motif=1 then flag_fru_motif=1;
    else flag_fru_motif=0;
    run;

proc freq data=motif;
    tables flag_long*flag_male_a_ind*flag_fru_a_motif/cmh;
    tables flag_long*flag_male_b_ind*flag_fru_b_motif/cmh;
    tables flag_long*flag_male_c_ind*flag_fru_c_motif/cmh;

    tables flag_long*flag_female_a_ind*flag_fru_a_motif/cmh;
    tables flag_long*flag_female_b_ind*flag_fru_b_motif/cmh;
    tables flag_long*flag_female_c_ind*flag_fru_c_motif/cmh;

    tables flag_long*flag_null_ind*flag_fru_a_motif/cmh;
    tables flag_long*flag_null_ind*flag_fru_b_motif/cmh;
    tables flag_long*flag_null_ind*flag_fru_c_motif/cmh;

    tables flag_long*flag_null_ind*flag_fru_motif/cmh;

    run;

proc freq data=motif;
    where flag_fru_a_motif=0 and flag_fru_c_motif=0;
    tables flag_long*flag_male_b_ind*flag_fru_b_motif/cmh;

    tables flag_long*flag_female_b_ind*flag_fru_b_motif/cmh;

    tables flag_long*flag_null_ind*flag_fru_b_motif/cmh;
    run;

proc freq data=motif;
    tables flag_fru_a_motif*(flag_fru_b_motif flag_fru_c_motif)/chisq;
    run;
