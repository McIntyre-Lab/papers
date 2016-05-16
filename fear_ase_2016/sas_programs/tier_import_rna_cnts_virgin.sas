/*******************************************************************************
* Filename: tier_prep_r101.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Going to use RNA reads form R101 and the DNA will be importated
* from the read_simulation.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Import all virgin samples from design file */
    proc sort data=CEGS.ase_design_file;
        by line mating_status rep;
        run;

    data design_file;
        set CEGS.ase_design_file;
        where mating_status eq 'V';
        run;

/* Import sam-compare RNA Counts */
    %macro import_sam_compare(line, ms, rep);
        %let name = "!MCLAB/cegs_ase_paper/pipeline_output/ase_counts_fb551_updated_fusions/ase_counts_&line._&ms.&rep..csv";

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
                retain line mating_status rep fusion_id prop_RNA LINE_RNA_TOTAL TESTER_RNA_TOTAL ASE_RNA_TOTAL BOTH_RNA_TOTAL flag_simulate;
                length line $7. ;
                set ase_counts;
                line = "&line";
                mating_status = "&ms";
                rep = "&rep";

                TOTAL_READS_COUNTED = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL,
                                          SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT, SAM_A_EXACT_SAM_B_INEXACT,
                                          SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                          SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER); 

                BOTH_RNA_TOTAL = sum(BOTH_EXACT, BOTH_INEXACT_EQUAL);

                LINE_RNA_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                 SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER);

                TESTER_RNA_TOTAL = sum(SAM_B_ONLY_EXACT, SAM_B_ONLY_SINGLE_INEXACT,
                                   SAM_B_EXACT_SAM_A_INEXACT, SAM_B_INEXACT_BETTER);

                ASE_RNA_TOTAL = sum(SAM_A_ONLY_EXACT, SAM_A_ONLY_SINGLE_INEXACT,
                                SAM_A_EXACT_SAM_B_INEXACT, SAM_A_INEXACT_BETTER, SAM_B_ONLY_EXACT,
                                SAM_B_ONLY_SINGLE_INEXACT, SAM_B_EXACT_SAM_A_INEXACT,
                                SAM_B_INEXACT_BETTER);

                if ASE_RNA_TOTAL > 0 then do;
                    prop_rna = LINE_RNA_TOTAL/ASE_RNA_TOTAL ;
                    flag_simulate = 1;
                end;
                else do;
                    prop_rna = 0.5;
                    flag_simulate = 0;
                end;

                keep line mating_status rep fusion_id prop_rna LINE_RNA_TOTAL TESTER_RNA_TOTAL ASE_RNA_TOTAL BOTH_RNA_TOTAL flag_simulate;
                run;

            proc append base=rna_cnts data=ase_counts2;
            run;

        %end;
    %mend;
    %iterdataset(dataset=design_file,function=%nrstr(%import_sam_compare(&line,&mating_status,&rep);));

/* Clean Up */
    proc datasets nolist;
        delete ase_counts;
        delete ase_counts2;
        delete design_file;
        run; quit;
