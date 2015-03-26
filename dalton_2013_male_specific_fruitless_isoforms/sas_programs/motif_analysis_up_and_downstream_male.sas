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

    %macro response_type(type);

        proc freq data=ind_rep_motif noprint;
            tables flag_male_&letter._&type.*flag_fru_&letter._motif /chisq expected outexpect out=fru_&letter._&type._motif;
            exact fisher chisq;
            output out=fru_&letter._&type._motif_test chisq fisher;
            run;

        data fru_&letter._&type._motif_test;
            set fru_&letter._&type._motif_test;
            length name $9.;
            name = "&type.";
            label XP2_FISH = 'Fisher raw P-value (2-tail)';
            keep name N _PCHI_ DF_PCHI XP_PCHI XP2_FISH;
            run;

    %mend response_type;
    %response_type(ind);
    %response_type(rep);

    data fru_&letter._motif_tests;
        retain sample name;
        length sample $5.;
        set Fru_&letter._ind_motif_test Fru_&letter._rep_motif_test;
        sample = "fru_&letter.";
        run;

    proc datasets nolist;
        delete fru_&letter._ind_motif_test;
        delete fru_&letter._rep_motif_test;
        run;

%mend enrichment_analysis;
%enrichment_analysis(a);
%enrichment_analysis(b);
%enrichment_analysis(c);

data FRU.fru_motif_test_up_down_male;
    set fru_a_motif_tests fru_b_motif_tests fru_c_motif_tests ;
    run;

/** Export Contengency Tables **/
%macro export_analysis(letter,outdir);

    %macro response_type(type);

        proc export data=fru_&letter._&type._motif 
                    outfile="&outdir./fru_&letter._&type._motif_table_up_and_down_male.csv" 
                    dbms=CSV replace label;
            putnames=yes;
            run;

    %mend response_type;
    %response_type(ind);
    %response_type(rep);

%mend export_analysis;
%export_analysis(a,!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);
%export_analysis(b,!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);
%export_analysis(c,!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);


/* Export test Statistics */
proc export data=fru.fru_motif_test_up_down_male 
            outfile="!MCLAB/arbeitman/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_motif_tests_up_and_down_male.csv" 
            dbms=CSV replace label;
    putnames=yes;
    run;

proc datasets nolist;
    delete FLAG_IND_REP FLAG_MOTIF FRU_A_IND_MOTIF FRU_A_MOTIF_TESTS
    FRU_A_REP_MOTIF FRU_B_IND_MOTIF FRU_B_MOTIF_TESTS FRU_B_REP_MOTIF
    FRU_C_IND_MOTIF FRU_C_MOTIF_TESTS FRU_C_REP_MOTIF IND_REP_MOTIF ;
    run;

