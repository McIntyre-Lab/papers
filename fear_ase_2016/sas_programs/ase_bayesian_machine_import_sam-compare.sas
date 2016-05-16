/********************************************************************************
* Import Sam Compare and build stack dataset
********************************************************************************/
/*
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname cegsqc '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/qc/sas_data';
libname design '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/design_files/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Import MACRO */
    %macro import_sam_compare(folder, line, mating_status, rep);
        %let name = "&folder./ase_counts_&line._w1118&line._&mating_status._&rep..csv";

        %local rc fileref;
        %let rc = %sysfunc(filename(fileref,&name));

        %if %sysfunc(fexist(&fileref)) %then %do;
            data WORK.ASE_COUNTS    ;
                infile &name delimiter = ',' MISSOVER DSD lrecl=32767  firstobs=2 ;
                informat FUSION_ID $9. ;
                informat BOTH_EXACT best32. ;
                informat BOTH_INEXACT_EQUAL best32. ;
                informat SAM_A_ONLY_EXACT best32. ;
                informat SAM_B_ONLY_EXACT best32. ;
                informat SAM_A_EXACT_SAM_B_INEXACT best32. ;
                informat SAM_B_EXACT_SAM_A_INEXACT best32. ;
                informat SAM_A_ONLY_SINGLE_INEXACT best32. ;
                informat SAM_B_ONLY_SINGLE_INEXACT best32. ;
                informat SAM_A_INEXACT_BETTER best32. ;
                informat SAM_B_INEXACT_BETTER best32. ;
                format FUSION_ID $9. ;
                format BOTH_EXACT best12. ;
                format BOTH_INEXACT_EQUAL best12. ;
                format SAM_A_ONLY_EXACT best12. ;
                format SAM_B_ONLY_EXACT best12. ;
                format SAM_A_EXACT_SAM_B_INEXACT best12. ;
                format SAM_B_EXACT_SAM_A_INEXACT best12. ;
                format SAM_A_ONLY_SINGLE_INEXACT best12. ;
                format SAM_B_ONLY_SINGLE_INEXACT best12. ;
                format SAM_A_INEXACT_BETTER best12. ;
                format SAM_B_INEXACT_BETTER best12. ;
                input
                    FUSION_ID $
                    BOTH_EXACT
                    BOTH_INEXACT_EQUAL
                    SAM_A_ONLY_EXACT
                    SAM_B_ONLY_EXACT
                    SAM_A_EXACT_SAM_B_INEXACT
                    SAM_B_EXACT_SAM_A_INEXACT
                    SAM_A_ONLY_SINGLE_INEXACT
                    SAM_B_ONLY_SINGLE_INEXACT
                    SAM_A_INEXACT_BETTER
                    SAM_B_INEXACT_BETTER
                    ;
                run;

            data ase_counts2;
                retain line mating_status rep fusion_id TOTAL_READS_COUNTED BOTH_TOTAL LINE_TOTAL TESTER_TOTAL ASE_TOTAL;
                length line $13. ;
                length mating_status $1. ;
                length rep $5. ;
                set ase_counts;
                line = "&line";
                mating_status = "&mating_status";
                rep = "&rep";

                TOTAL_READS_COUNTED = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL,
                                          SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT, SAM_A_EXACT_SAM_B_INEXACT,
                                          SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                          SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER); 

                BOTH_TOTAL = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL);

                LINE_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                 SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER);

                TESTER_TOTAL = sum(SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                   SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER);

                ASE_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT,
                                SAM_B_ONLY_SINGLE_INEXACT, SAM_B_EXACT_SAM_A_INEXACT,
                                SAM_B_INEXACT_BETTER);

                drop BOTH_EXACT BOTH_INEXACT_EQUAL SAM_A_ONLY_EXACT SAM_B_ONLY_EXACT SAM_A_EXACT_SAM_B_INEXACT SAM_B_EXACT_SAM_A_INEXACT
                SAM_A_ONLY_SINGLE_INEXACT SAM_B_ONLY_SINGLE_INEXACT SAM_A_INEXACT_BETTER SAM_B_INEXACT_BETTER;

                run;

            proc append base=all_ase data=ase_counts2;
            run;
        %end;
    %mend;

/* Import ase counts for complete folder */
    %let folder=!MCLAB/cegs_ase_paper/pipeline_output/ase_counts_fb551_updated_fusions;

    proc sort data=CEGS.complete_design_by_rep;
        by line mating_status rep;
        run;

    proc sort data=CEGS.flag_contaminated;
        by line mating_status rep;
        run;

    * Do not import contaminated samples;
    data design bad;
        merge CEGS.complete_design_by_rep CEGS.flag_contaminated;
        by line mating_status rep;
        if flag_contaminated = 0 then output design;
        else output bad;
        drop flag_contaminated;
        run;

    %iterdataset(dataset=design,function=%nrstr(%import_sam_compare(&folder, &line, &mating_status, &rep);));


/* Make Permanent dataset */
    data CEGS.ase_counts_for_bayesian;
        set all_ase;
        run;

/* Clean Up */
    proc datasets nolist;
        delete ase_counts;
        delete ase_counts2;
        delete design;
        delete bad;
        delete all_ase;
        run; quit;
