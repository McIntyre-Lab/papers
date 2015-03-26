/********************************************************************************
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';
* libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Motif Enrichment Analysis */

data flag_ind_rep;
    set FRU.flag_ind_rep;
    run;

data flag_motif;
    set FRU.motif_flags_and_cnts;
    if flag_fru_a_motif = 1 or flag_flag_fru_b_motif = 1 or flag_fru_c_motif = 1 then flag_motif = 1;
    else flag_motif = 0;
    keep primary_fbgn flag_fru_a_motif flag_fru_b_motif flag_fru_c_motif flag_motif;
    run;

proc sort data=flag_ind_rep;
    by primary_fbgn;
    run;

proc sort data=flag_motif;
    by primary_fbgn;
    run;

data ind_rep_motif;
    merge flag_ind_rep flag_motif;
    by primary_fbgn;
    run;


%macro enrichment_analysis(letter);

        proc freq data=ind_rep_motif noprint;
            tables flag_male_&letter.*flag_fru_&letter._motif /chisq expected outexpect out=fru_male_&letter._motif;
            exact fisher chisq;
            output out=fru_&letter._motif_test chisq fisher;
            run;

        data fru_&letter._motif_test;
            length sample $5.;
            retain sample;
            set fru_&letter._motif_test;
            label XP2_FISH = 'Fisher raw P-value (2-tail)';
            sample = "fru_&letter";
            keep sample N _PCHI_ DF_PCHI XP_PCHI XP2_FISH;
            run;

%mend enrichment_analysis;
%enrichment_analysis(a);
%enrichment_analysis(b);
%enrichment_analysis(c);

data FRU.fru_motif_test_up_down_male_v2;
    set fru_a_motif_test fru_b_motif_test fru_c_motif_test ;
    run;

/** Export Contengency Tables **/
%macro export_analysis(letter,outdir);

        proc export data=fru_male_&letter._motif 
                    outfile="&outdir./fru_&letter._motif_table_up_and_down_male_v2.csv" 
                    dbms=CSV replace label;
            putnames=yes;
            run;

%mend export_analysis;
%export_analysis(a,!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);
%export_analysis(b,!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);
%export_analysis(c,!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);


/* Export test Statistics */
proc export data=FRU.fru_motif_test_up_down_male_v2
            outfile="!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_male_v2.csv" 
            dbms=CSV replace label;
    putnames=yes;
    run;

proc datasets nolist;
    delete FLAG_IND_REP FLAG_MOTIF FRU_A_MOTIF_TEST FRU_B_MOTIF_TEST
    FRU_C_MOTIF_TEST FRU_MALE_A_MOTIF FRU_MALE_B_MOTIF FRU_MALE_C_MOTIF
    IND_REP_MOTIF;
    run;
