/*******************************************************************************
* Filename: tier_merge_export_r101.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Combine the RNA and simulated DNA counts for r101 and export
* them for the TIER simulation.
*
*******************************************************************************/

/* Libraries
libname CEGS '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/* Grab fusions that have ASE reads in Virgins */
    * lets make life easy and only use 3 reps;
    data filter;
        set rna_cnts;
        where flag_simulate eq 1 and mating_status eq 'V';
        if rep le 3;
        keep line fusion_id rep LINE_RNA_TOTAL TESTER_RNA_TOTAL;
        run;

    proc sort data=filter;
        by line fusion_id rep;
        run;

/* Make Wide Dataset */
    data wide;
        set filter;
        by line fusion_id rep;
        keep line fusion_id line_RNA_Rep: tester_RNA_Rep:;
        retain line fusion_id line_RNA_Rep1-line_RNA_Rep3 tester_RNA_Rep1-tester_RNA_Rep3;
        array linecnts(1:3) line_RNA_Rep1-line_RNA_Rep3;
        array testercnts(1:3) tester_RNA_Rep1-tester_RNA_Rep3;
        if first.fusion_id then do;
            do i = 1 to 3;
                linecnts(i) = .;
                testercnts(i) = .;
            end;
        end;
        linecnts(trim(rep)) = line_RNA_total;
        testercnts(trim(rep)) = tester_RNA_total;
        if last.fusion_id then output;
        run;

/* Drop fusions that don't have data for all reps or that have 0's */
    data filter2;
        set wide;
        if line_RNA_Rep1 = . or
        line_RNA_Rep2 = . or
        line_RNA_Rep3 = . or
        tester_RNA_Rep1 = . or
        tester_RNA_Rep2 = . or
        tester_RNA_Rep3 = . then delete;
        if (line_RNA_Rep1 = 0 and line_RNA_Rep2 = 0 and line_RNA_Rep3 = 0) or
        (tester_RNA_Rep1 = 0 and tester_RNA_Rep2 = 0 and tester_RNA_Rep3 = 0) then delete;
        run;

/* Merge DNA and RNA */
    proc sort data=filter2;
        by line fusion_id;
        run;

    proc sort data=dna_cnts;
        by line fusion_id;
        run;

    data merged;
        merge filter2 (in=in1) dna_cnts (in=in2 );
        by line fusion_id;
        if in1 and in2;
        drop ASE_DNA_TOTAL BOTH_DNA_TOTAL;
        if LINE_DNA_TOTAL le 10 then delete;
        run;

/* Export by line */
    data design_file;
        set CEGS.ase_design_file;
        keep line;
        run;

    proc sort data=design_file nodupkey;
        by line;
        run;

    %macro exportLine(line);
        data out;
            set merged;
            where line eq '&line';
            run;

        proc export data=out outfile='!MCLAB/cegs_ase_paper/pipeline_output/typeI_error/input/&line._RNA_sim_DNA_cnts.csv' dbms=csv replace;
            putnames=yes;
            run;
    %mend;

    %iterdataset(dataset=design_file, function=%nrstr(%exportLine(&line);));


/* Clean UP */
    proc datasets nolist;
        delete filter;
        delete filter2;
        delete merged;
        delete rna_cnts;
        delete wide;
        delete out;
        delete design_file;
        delete dna_cnts;
        run; quit;

