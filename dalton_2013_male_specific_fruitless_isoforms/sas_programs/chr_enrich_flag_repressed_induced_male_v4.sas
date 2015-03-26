/********************************************************************************
* This script creates a flag for genes that are repressed or induced by Fru.
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';
* libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/** Rename var1 to primary_fbgn for merging **/
%macro rename_var1(in_data);
data &in_data;
    set &in_data;
    rename var1=primary_fbgn;
    run;
%mend rename_var1;

%rename_var1(fru.male_induced_union);
%rename_var1(fru.male_repressed_union);


/** Sort by primary_fbgn for merging **/
%macro sort_fbgn(in_data);
proc sort data=&in_data;
    by primary_fbgn;
    run;
%mend sort_fbgn;

%sort_fbgn(fru.male_induced_union);
%sort_fbgn(fru.male_repressed_union);
%sort_fbgn(dmel530.symbol2fbgn);
%sort_fbgn(dmel530.fbgn2coord);

/** Merge Datasets and Flag Induced and Repressed **/

data some_merged;
    merge dmel530.symbol2fbgn (in=in1) dmel530.fbgn2coord (in=in2);
    by primary_fbgn;
    if in2;
    run;

%sort_fbgn(some_merged);

data all_merged oops oops2;
    merge fru.male_induced_union (in=in1) fru.male_repressed_union (in=in2) some_merged (in=in3);
    by primary_fbgn;
    if in1 then flag_induced = 1; else flag_induced = 0;
    if in2 then flag_repressed = 1; else flag_repressed = 0;
    if flag_a_ind = " " then flag_a_ind = 0;
    if flag_b_ind = " " then flag_b_ind = 0;
    if flag_c_ind = " " then flag_c_ind = 0;
    if flag_a_rep = " " then flag_a_rep = 0;
    if flag_b_rep = " " then flag_b_rep = 0;
    if flag_c_rep = " " then flag_c_rep = 0;
    if in1 and not in3 then output oops;
    else if in2 and not in3 then output oops;
    else output all_merged;
	
	if in3 and not (in2 or in1) then output oops2 ;
    run;   *oops has 0 ;
    * oops2 has 14179, There are 5 genes that are induced in one and repressed in another ;
    * so 267 + 462 = 729 - 5 = 724 ;
    * 14903 - 724 = 14179;

/* Checks

    * Check overlap of induced and repressed;
    proc freq data=all_merged;
        tables flag_induced*flag_repressed;
        run; * There are 5 obs that are both induced and repressed;

    * Check that all my chromosomes are present;
    proc freq data=all_merged;
        tables chrom;
        run;

*/

data fru.flag_x_induced_repressed_male;
    retain primary_fbgn symbol chrom start end strand;
    set all_merged;
    if chrom eq "X" or chrom eq "XHet" then flag_x_chrom = 1; else flag_x_chrom = 0;
    if chrom eq "2L" or chrom eq "2LHet" then flag_2L_chrom = 1; else flag_2L_chrom = 0;
    if chrom eq "2R" or chrom eq "2RHet" then flag_2R_chrom = 1; else flag_2R_chrom = 0;
    if chrom eq "3L" or chrom eq "3LHet" then flag_3L_chrom = 1; else flag_3L_chrom = 0;
    if chrom eq "3R" or chrom eq "3RHet" then flag_3R_chrom = 1; else flag_3R_chrom = 0;
    if chrom eq "4" then flag_4_chrom = 1; else flag_4_chrom = 0;
    run; * 14903 obs ;


data fru.flag_ind_rep_w_het_male;
    retain primary_fbgn symbol chrom start end strand;
    set all_merged;
    if chrom eq "X" then flag_x_chrom = 1; else flag_x_chrom = 0;
    if chrom eq "XHet" then flag_xhet_chrom = 1; else flag_xhet_chrom = 0;
    if chrom eq "2L" then flag_2L_chrom = 1; else flag_2L_chrom = 0;
    if chrom eq "2LHet" then flag_2Lhet_chrom = 1; else flag_2Lhet_chrom = 0;
    if chrom eq "2R" then flag_2R_chrom = 1; else flag_2R_chrom = 0;
    if chrom eq "2RHet" then flag_2Rhet_chrom = 1; else flag_2Rhet_chrom = 0;
    if chrom eq "3L" then flag_3L_chrom = 1; else flag_3L_chrom = 0;
    if chrom eq "3LHet" then flag_3Lhet_chrom = 1; else flag_3Lhet_chrom = 0;
    if chrom eq "3R" then flag_3R_chrom = 1; else flag_3R_chrom = 0;
    if chrom eq "3RHet" then flag_3Rhet_chrom = 1; else flag_3Rhet_chrom = 0;
    if chrom eq "4" then flag_4_chrom = 1; else flag_4_chrom = 0;
    run; * 14903 obs ;
