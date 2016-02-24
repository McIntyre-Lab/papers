/*******************************************************************************
* Filename: Import_indel_vcf_files_v2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Import split INDEL vcf files into SAS
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase557/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/*
data design_file;
    set CEGS.genotype_list;
    run;
*/

*Import the VCF files for each line;
%macro Import_vcf (ID, first);
    data WORK.&ID._indel_vcf;
        infile "!HOME/sandbox/cegs_ase_paper/indels_by_sample/&ID..vcf" delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=&first ;
        informat CHROM $32. ;
        informat POS best32. ;
        informat ID_&ID. best32. ;
        informat REF_&ID. $300. ;
        informat ALT_&ID. $300. ;
        informat QUAL_&ID. best32. ;
        informat FILTER_&ID. $36. ;
        informat INFO_&ID. $262. ;
        informat FORMAT_&ID. $18. ;
        informat &ID. $300. ;
        format CHROM $32. ;
        format POS best32. ;
        format ID_&ID. best32. ;
        format REF_&ID. $300. ;
        format ALT_&ID. $300. ;
        format QUAL_&ID. best32. ;
        format FILTER_&ID. $36. ;
        format INFO_&ID. $262. ;
        format FORMAT_&ID. $18. ;
        format &ID. $300. ;
        input
            CHROM $
            POS
            ID_&ID.
            REF_&ID. $
            ALT_&ID. $
            QUAL_&ID.
            FILTER_&ID. $
            INFO_&ID. $
            FORMAT_&ID. $
            &ID. $
        ;
        run;
%mend;

%Import_vcf(w1118, 61);
%iterdataset(dataset=design_file, function=%nrstr(%Import_vcf(&line, 61);));

