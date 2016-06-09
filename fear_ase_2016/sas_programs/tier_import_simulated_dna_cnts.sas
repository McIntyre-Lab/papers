/*******************************************************************************
* Filename: tier_import_simulated_dna_cnts_r101.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Going to use DNA counts form the R101 read simulation for simulating data
* for doing TIER
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Make list of lines */
    data design_file;
        set CEGS.ase_design_file;
        keep line;
        run;

    proc sort data=design_file nodupkey;
        by line;
        run;

/* Import sam-compare DNA Counts */
    %macro import_sam_compare(line);
        %let name = "!MCLAB/cegs_ase_paper/pipeline_output/qsim_bayesian/ase_counts_fb551/ase_counts_&line..csv";

        %local rc fileref;
        %let rc = %sysfunc(filename(fileref,&name));
        %if %sysfunc(fexist(&fileref)) %then %do;
            data WORK.ASE_COUNTS    ;
                infile &name delimiter = ',' MISSOVER DSD lrecl=32767  firstobs=2 ;
                informat FUSION_ID $20. ;
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
                format FUSION_ID $20. ;
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
                retain line fusion_id prop_DNA LINE_DNA_TOTAL TESTER_DNA_TOTAL ASE_DNA_TOTAL BOTH_DNA_TOTAL ;
                length line $7. ;
                set ase_counts;
                line = "&line";

                TOTAL_READS_COUNTED = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL,
                                          SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT, SAM_A_EXACT_SAM_B_INEXACT,
                                          SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                          SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER); 

                BOTH_DNA_TOTAL = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL);

                LINE_DNA_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                 SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER);

                TESTER_DNA_TOTAL = sum(SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                   SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER);

                ASE_DNA_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT,
                                SAM_B_ONLY_SINGLE_INEXACT, SAM_B_EXACT_SAM_A_INEXACT,
                                SAM_B_INEXACT_BETTER);

                if ASE_DNA_TOTAL > 0 then do;
                    prop_dna = LINE_DNA_TOTAL/ASE_DNA_TOTAL ;
                end;
                else do;
                    prop_dna = 0.5;
                end;

                keep line fusion_id prop_dna LINE_DNA_TOTAL TESTER_DNA_TOTAL ASE_DNA_TOTAL BOTH_DNA_TOTAL ;
                run;

            proc append base=dna_cnts data=ase_counts2;
            run;

        %end;
    %mend;
    %iterdataset(dataset=design_file, function=%nrstr(%import_sam_compare(&line);));

/* Clean Up */
    proc datasets nolist;
        delete ASE_COUNTS;
        delete ASE_COUNTS2;
        delete design_file;
        run; quit;
