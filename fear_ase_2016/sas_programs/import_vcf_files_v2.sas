/*******************************************************************************
* Filename: Import_vcf_files_v2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Import split vcf files into SAS
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase557/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/*
data design_file;
    set CEGS.genotype_list;
    run;
*/

*Import the snp VCF files for each line;
%macro import_snp_vcf(ID, first);
    data WORK.&ID._vcf ;
        infile "!HOME/sandbox/cegs_ase_paper/snps_by_sample/&ID..vcf" delimiter='09'x MISSOVER DSD lrecl=32767 firstobs=&first ;
        informat CHROM $32. ;
        informat POS best32. ;
        informat ID_&ID. best32. ;
        informat REF_&ID. $300. ;
        informat ALT_&ID. $300. ;
        informat QUAL_&ID. best32. ;
        informat FILTER_&ID. $36. ;
        informat INFO_&ID. $262. ;
        informat FORMAT_&ID. $18. ;
        informat &ID. $74. ;
        format CHROM $32. ;
        format POS best32. ;
        format ID_&ID. best32. ;
        format REF_&ID. $300. ;
        format ALT_&ID. $300. ;
        format QUAL_&ID. best32. ;
        format FILTER_&ID. $36. ;
        format INFO_&ID. $262. ;
        format FORMAT_&ID. $18. ;
        format &ID. $74. ;
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

%import_snp_vcf(w1118, 58);
%iterdataset(dataset=design_file, function=%nrstr(%import_snp_vcf(&line, 58);));

/* Check number of obs */
%macro count_obs (ID) ;
    title &id;
    proc contents data = WORK.&ID._vcf ; 
    run;
    title;
%mend ;
/*
%count_obs(w1118);
%iterdataset(dataset=design_file, function=%nrstr(%count_obs(&line);));
*/

/*
proc freq data = WORK.r101_vcf ;
tables ref_r101*alt_r101;
run;

proc freq data = WORK.r109_vcf ;
tables ref_r109*alt_r109;
run;
*/
