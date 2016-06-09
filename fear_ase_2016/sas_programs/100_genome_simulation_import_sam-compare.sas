/*******************************************************************************
* Filename: 100_genome_simulation_import_sam-compare.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Import sam-compare results
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Import MACRO */
    %macro import_sam_compare();
        %do ii = 1 %to 100;
            %let name = "!MCLAB/cegs_ase_paper/pipeline_output/100_genome_simulation/ase_counts_fb551/ase_counts_Line&ii..csv";

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
                    retain line mating_status rep fusion_id TOTAL_READS_COUNTED BOTH_TOTAL LINE_TOTAL TESTER_TOTAL ASE_TOTAL;
                    length line $7. ;
                    set ase_counts;
                    line = "Line&ii";

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

                    if ASE_TOTAL > 0 then PROP_TESTER = TESTER_TOTAL / ASE_TOTAL;
                    else PROP_TESTER = .;

                    drop BOTH_EXACT BOTH_INEXACT_EQUAL SAM_A_ONLY_EXACT SAM_B_ONLY_EXACT SAM_A_EXACT_SAM_B_INEXACT SAM_B_EXACT_SAM_A_INEXACT
                    SAM_A_ONLY_SINGLE_INEXACT SAM_B_ONLY_SINGLE_INEXACT SAM_A_INEXACT_BETTER SAM_B_INEXACT_BETTER;

                    run;

                proc append base=all_ase data=ase_counts2;
                run;
            %end;
        %end;
    %mend;

/* Import ase counts */
    %import_sam_compare();
    proc datasets nolist;
        delete ase_counts;
        delete ase_counts2;
        run; quit;
