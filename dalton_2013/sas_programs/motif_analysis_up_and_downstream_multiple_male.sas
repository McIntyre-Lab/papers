/********************************************************************************
********************************************************************************/

* libname fru '!MCLAB/arbeitman_fru_network/sasdata';
* libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

/* Multiple Motif Enrichment Analysis */
%macro enrichment_analysis(letter);

    %macro response_type(type);

        proc freq data=fru_&letter._motif noprint;
            tables flag_&letter._&type.*flag_multi_motif /chisq expected outexpect out=fru_&letter._&type._multi_motif;
            exact fisher chisq;
            output out=fru_&letter._&type._multi_motif_test chisq fisher;
            run;

        data fru_&letter._&type._multi_motif_test;
            set fru_&letter._&type._multi_motif_test;
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
        set Fru_&letter._ind_multi_motif_test Fru_&letter._rep_multi_motif_test;
        sample = "fru_&letter.";
        run;

    data fru_&letter._multi_motif_tests;
        retain sample name;
        length sample $5.;
        set Fru_&letter._ind_multi_motif_test Fru_&letter._rep_multi_motif_test;
        sample = "fru_&letter.";
        run;

    proc datasets nolist;
        delete fru_&letter._ind_motif_test;
        delete fru_&letter._rep_motif_test;
        delete fru_&letter._ind_multi_motif_test;
        delete fru_&letter._rep_multi_motif_test;
        run;

%mend enrichment_analysis;
%enrichment_analysis(a);
%enrichment_analysis(b);
%enrichment_analysis(c);

data fru.fru_multi_test_up_down_male;
    set fru_a_multi_motif_tests fru_b_multi_motif_tests fru_c_multi_motif_tests ;
    run;

/** Export Contengency Tables **/
%macro export_analysis(letter,outdir);

    %macro response_type(type);

        proc export data=fru_&letter._&type._multi_motif 
                    outfile="&outdir./fru_&letter._&type._multi_motif_table_up_and_down_male.csv" 
                    dbms=CSV replace label;
            putnames=yes;
            run;

    %mend response_type;
    %response_type(ind);
    %response_type(rep);

%mend export_analysis;
%export_analysis(a,!MCLAB/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);
%export_analysis(b,!MCLAB/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);
%export_analysis(c,!MCLAB/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down);


/* Export test Statistics */
proc export data=fru.fru_multi_test_up_down_male 
            outfile="!MCLAB/arbeitman_fru_network/reports_internal/motif_analysis/motif_results_up_and_down/fru_multi_motif_tests_up_and_down_male.csv" 
            dbms=CSV replace label;
    putnames=yes;
    run;

