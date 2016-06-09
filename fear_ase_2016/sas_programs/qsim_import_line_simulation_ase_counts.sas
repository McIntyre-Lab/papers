/*******************************************************************************
* Filename: qsim_import_line_simulation_ase_counts.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: This script does the following: (1) Imports sam-compare results
* from simulation from lines (2) calcualtes qsim_line = line_total / ASE_total
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Import MACRO */
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
                retain line fusion_id qsim_line LINE_SIM_TOTAL TESTER_SIM_TOTAL ASE_TOTAL BOTH_TOTAL;
                length line $5. ;
                set ase_counts;
                line = "&line";

                TOTAL_READS_COUNTED = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL,
                                          SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT, SAM_A_EXACT_SAM_B_INEXACT,
                                          SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                          SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER); 

                BOTH_TOTAL = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL);

                LINE_SIM_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                 SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER);

                TESTER_SIM_TOTAL = sum(SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                   SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER);

                ASE_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT,
                                SAM_B_ONLY_SINGLE_INEXACT, SAM_B_EXACT_SAM_A_INEXACT,
                                SAM_B_INEXACT_BETTER);

                if ASE_TOTAL > 0 then qsim_tester = TESTER_SIM_TOTAL/ASE_TOTAL ;
                else qsim_line = 0.5;

                keep line fusion_id qsim_tester LINE_SIM_TOTAL TESTER_SIM_TOTAL ASE_TOTAL BOTH_TOTAL;

                run;

            proc append base=all_ase data=ase_counts2;
            run;
        %end;
    %mend;
    %iterdataset(dataset=design_file,function=%nrstr(%import_sam_compare(&line);));

    proc sort data=all_ase;
        by fusion_id;
        run;

    data CEGS.ase_qsim_tester;
        set all_ase;
        keep line fusion_id qsim_tester;
        run;

    proc datasets nolist;
        delete ase_counts;
        delete ase_counts2;
        delete all_ase;
        run; quit;
