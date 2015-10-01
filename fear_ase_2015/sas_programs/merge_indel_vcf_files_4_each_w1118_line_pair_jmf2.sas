/*******************************************************************************
* Filename: Merge_vcf_files_4_each_w1118_line_pair_jmf2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper';
libname dmel '!MCLAB/useful_dmel_data/flybase557/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/*
data design_file;
    set CEGS.genotype_list;
    run;
*/

%macro merge_w11182line_vcf (ID);

    proc sort data = WORK.w1118_indel_vcf;
        by chrom pos;
        run;

    proc sort data = WORK.&ID._indel_vcf;
        by chrom pos;
        run;

    data w1118_2_&ID._indel_vcf;
        merge WORK.w1118_indel_vcf (in=in1) WORK.&ID._indel_vcf (in=in2);
        by chrom pos;
        if in1 and in2 then VCF_in_both = 1;
        else VCF_in_both = 0;
        if in1 and not in2 then VCF_in_w1118_only = 1;
        else VCF_in_w1118_only = 0;
        if not in1 and in2 then VCF_in_&ID._only = 1;
        else VCF_in_&ID._only = 0;
        run;

    *Flag Hets;
    data w1118_2_&ID._indel_vcf;
        set w1118_2_&ID._indel_vcf;
        find_w1118_het = count(w1118,'Het');
        if find_w1118_het >= 1 then w1118_ishet = 1;
        else w1118_ishet = 0;
        find_&ID._het = count(&ID.,'Het');
        if find_&ID._het >= 1 then &ID._ishet = 1;
        else &ID._ishet = 0;
        drop find_w1118_het find_&ID._het;
        run;

    /* Clean up */
    proc datasets nolist;
        delete &ID._indel_vcf;
        run; quit;
%mend;
%iterdataset(dataset=design_file, function=%nrstr(%merge_w11182line_vcf(&line);));

proc datasets nolist;
    delete w1118_indel_vcf;
    run; quit;
