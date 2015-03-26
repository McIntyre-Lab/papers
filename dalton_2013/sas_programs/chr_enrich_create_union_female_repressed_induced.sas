/********************************************************************************
* I need to create the union between Fru mutants for genes that are repressed
* or induced.
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';

/** Remove Duplicates from gene lists  **/
%macro remove_dups(in_data);
proc sort data=&in_data(where=(var1 ^? ';')) nodups out=&in_data._nodups;
    by var1;
    run;
%mend remove_dups;

%remove_dups(FruA_ind); * 111 obs;
%remove_dups(FruA_rep); * 183 obs;
%remove_dups(FruB_ind); * 117 obs;
%remove_dups(FruB_rep); * 237 obs;
%remove_dups(FruC_ind); * 167 obs;
%remove_dups(FruC_rep); * 198 obs;

/** Create Union of Induced Genes  **/

data fru.female_induced_union;
    merge FruA_ind_nodups (in=in1) FruB_ind_nodups (in=in2) FruC_ind_nodups (in=in3);
    by var1;
    if in1 then flag_a_ind = 1; else flag_a_ind = 0;
    if in2 then flag_b_ind = 1; else flag_b_ind = 0;
    if in3 then flag_c_ind = 1; else flag_c_ind = 0;
    run; * 267 obs;

/** Create Union of Repressed Genes  **/

data fru.female_repressed_union;
    merge FruA_rep_nodups (in=in1) FruB_rep_nodups (in=in2) FruC_rep_nodups (in=in3);
    by var1;
    if in1 then flag_a_rep = 1; else flag_a_rep = 0;
    if in2 then flag_b_rep = 1; else flag_b_rep = 0;
    if in3 then flag_c_rep = 1; else flag_c_rep = 0;
    run; * 462 obs;

/*
    proc sort data=fru.female_induced_union nodups out=test;
        by var1;
        run; *no dups;

    proc sort data=fru.female_repressed_union nodups out=test;
        by var1;
        run; *no dups;

*/
