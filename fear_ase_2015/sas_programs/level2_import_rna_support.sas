/*******************************************************************************
* Filename: level2_import_rna_support.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Import the summarized RNA counts for each genotype
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

%macro merge_RNA_cnts(ID);
    %* Import RNA counts;
    data WORK.&ID._COUNTS    ;
        infile "!HOME/sandbox/cegs_ase_paper/ase_masked_aln_rna_support/&ID..csv" delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat chrom $32. ;
        informat pos best32. ;
        informat ref $300. ;
        informat alt $300. ;
        informat dnaCnt best32. ;
        informat altRnaCnt best32. ;
        informat refRnaCnt best32. ;
        informat totRnaCnt best32. ;
        format chrom $32. ;
        format pos best32. ;
        format ref $300. ;
        format alt $300. ;
        format dnaCnt best12. ;
        format altRnaCnt best12. ;
        format refRnaCnt best12. ;
        format totRnaCnt best12. ;
        input
            chrom $
            pos
            ref $
            alt $
            dnaCnt
            altRnaCnt
            refRnaCnt
            totRnaCnt
            ;
        run;

%mend;
%iterdataset(dataset=design_file, function=%nrstr(%merge_RNA_cnts(&line);));
