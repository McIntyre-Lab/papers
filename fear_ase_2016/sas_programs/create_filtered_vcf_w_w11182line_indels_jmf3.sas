/*******************************************************************************
* Filename: Create_filtered_vcf_w_w11182line_SNPs_jmf2.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper/indels';
libname dmel '!MCLAB/useful_dmel_data/flybase557/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/********************************************************************************/
/* Categorize SNP classes                                                       */
/* Flag SNPS for Level1 outputting                                              */
/* Output VCF files for creating masked alignments                              */
/* do not delete THUMP.w1118_2_&ID._vcf_flags2 !!                              */
/********************************************************************************/

/*
data design_file;
    set CEGS.genotype_list;
    run;
*/

%macro flag_w11182line_vcf (ID);
    data THUMP.w1118_2_&ID._indel_vcf_flags2;
        set WORK.w1118_2_&ID._indel_vcf;

        *Check if the reference base is the same;
        if VCF_in_both = 1 and REF_w1118 ~= REF_&ID then REF_NOT_IDENTICAL = 1;
        else REF_NOT_IDENTICAL = 0;
        if VCF_in_both ~= 1 then REF_NOT_IDENTICAL = ".";

        *Case 1, either position is heterozygous;
        if w1118_ishet = 1 or &ID._ishet = 1 then do;
            Output_Lvl1_MSKD = 0;
            Output_Lvl1_UPD = 0;
            case = 1;
        end;

        *Case 2, both ATL calls are the reference;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 = ALT_&ID and ALT_w1118 = "" and ALT_&ID = "" then do;
            Output_Lvl1_MSKD = 0;
            Output_Lvl1_UPD = 0;
            case = 2;
        end;

        *Case 3, both ATL calls are the same;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 = ALT_&ID and ALT_w1118 ~= "" and ALT_&ID ~= "" then do;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 3;
        end;

        *Case 4, only ALT_w1118 or ALT_&ID is the reference;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 ~= ALT_&ID and ALT_w1118 = "" then do;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 4;
        end;

        if w1118_ishet = 0 and VCF_in_both = 1 and ALT_w1118 ~= ALT_&ID and ALT_&ID = "" then do;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 4;
        end;

        *Case 5, both ALT calls are different from reference and different from each other ;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 ~= ALT_&ID and ALT_&ID ~= "" and ALT_w1118 ~= "" then do;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 5;
        end;

        *Case 6, the ALT call is only in one genotype, is different than the reference, and is not a het;
        if &ID._ishet = 0 and VCF_in_&ID._only = 1 and ALT_&ID ~= "" then do;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 6;
        end;

        if w1118_ishet = 0 and VCF_in_w1118_only = 1 and ALT_w1118 ~= "" then do;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 6;
        end;

        *Case 7, the ALT call is only in one genotype, is the same as the reference, and is not a het;
        if VCF_in_&ID._only  = 1 and ALT_&ID = "" then do;
            Output_Lvl1_MSKD = 0;
            Output_Lvl1_UPD = 0;
            case = 7;
        end;

        if VCF_in_w1118_only = 1 and ALT_w1118 = "" then do;
            Output_Lvl1_MSKD = 0;
            Output_Lvl1_UPD = 0;
            case = 7;
        end;
        run;

        proc datasets nolist;
            delete w1118_2_&ID._indel_vcf;
            run;quit;

%mend;
%iterdataset(dataset=design_file, function=%nrstr(%flag_w11182line_vcf(&line);));


*Output a VCF file for masking ;
%macro output_vcf(ID) ;
    data out;
        set THUMP.w1118_2_&ID._indel_vcf_flags2;
        where output_lvl1_MSKD=1;
        run;

    data w1118;
        retain chrom pos MSKD_ID MSKD_REF MSKD_ALT MSKD_QUAL MSKD_FILTER MSKD_INFO MSKD_FORMAT MSKD_GENOTYPE;
        set out;
        if ALT_w1118 ne '';
        MSKD_ID = ID_w1118;
        MSKD_REF = REF_w1118;
        MSKD_ALT = ALT_w1118; 
        MSKD_QUAL = QUAL_w1118;
        MSKD_FILTER = FILTER_w1118 ;
        MSKD_INFO = INFO_w1118;
        MSKD_FORMAT = FORMAT_w1118;
        MSKD_GENOTYPE = w1118;
        keep chrom pos MSKD_ID MSKD_REF MSKD_ALT MSKD_QUAL MSKD_FILTER MSKD_INFO MSKD_FORMAT MSKD_GENOTYPE;
        run;

    data _null_;
        file "!HOME/sandbox/cegs_ase_paper/ase_lvl1_filtered_vcf_files/w1118_w11182&ID._lvl1_indel.vcf"
        delimiter='09'x DSD DROPOVER lrecl=32767;
        if _n_ = 1 then do;
            put '#CHROM' '09'x 'POS' '09'x 'ID' '09'x 'REF' '09'x 'ALT' '09'x 'QUAL' '09'x 
                'FILTER' '09'x 'INFO' '09'x 'FORMAT' '09'x 'GENOTYPE';
        end;
        set  w1118   end=EFIEOD;
        format chrom $32. ;
        format pos best12. ;
        format MSKD_ID best12. ;
        format MSKD_REF $300. ;
        format MSKD_ALT $300. ;
        format MSKD_QUAL best32. ;
        format MSKD_FILTER $36. ;
        format MSKD_INFO $262. ;
        format MSKD_FORMAT $18. ;
        format MSKD_GENOTYPE $300. ;
        do;
            put chrom $ @;
            put pos @;
            put MSKD_ID $ @;
            put MSKD_REF $ @;
            put MSKD_ALT $ @;
            put MSKD_QUAL @;
            put MSKD_FILTER $ @;
            put MSKD_INFO $ @;
            put MSKD_FORMAT $ @;
            put MSKD_GENOTYPE $ ;
            ;
        end;
        run;

    data line;
        retain chrom pos MSKD_ID MSKD_REF MSKD_ALT MSKD_QUAL MSKD_FILTER MSKD_INFO MSKD_FORMAT MSKD_GENOTYPE;
        set out;
        if ALT_&ID. ne '';
        MSKD_ID = ID_&ID.;
        MSKD_REF = REF_&ID.;
        MSKD_ALT = ALT_&ID.; 
        MSKD_QUAL = QUAL_&ID.;
        MSKD_FILTER = FILTER_&ID. ;
        MSKD_INFO = INFO_&ID.;
        MSKD_FORMAT = FORMAT_&ID.;
        MSKD_GENOTYPE = &ID;
        keep chrom pos MSKD_ID MSKD_REF MSKD_ALT MSKD_QUAL MSKD_FILTER MSKD_INFO MSKD_FORMAT MSKD_GENOTYPE;
        run;

    data _null_;
        file "!HOME/sandbox/cegs_ase_paper/ase_lvl1_filtered_vcf_files/&ID._w11182&ID._lvl1_indel.vcf"
        delimiter='09'x DSD DROPOVER lrecl=32767;
        if _n_ = 1 then do;
            put '#CHROM' '09'x 'POS' '09'x 'ID' '09'x 'REF' '09'x 'ALT' '09'x 'QUAL' '09'x 
                'FILTER' '09'x 'INFO' '09'x 'FORMAT' '09'x 'GENOTYPE';
        end;
        set  line   end=EFIEOD;
        format chrom $32. ;
        format pos best12. ;
        format MSKD_ID best12. ;
        format MSKD_REF $300. ;
        format MSKD_ALT $300. ;
        format MSKD_QUAL best32. ;
        format MSKD_FILTER $36. ;
        format MSKD_INFO $262. ;
        format MSKD_FORMAT $18. ;
        format MSKD_GENOTYPE $300. ;
        do;
            put chrom $ @;
            put pos @;
            put MSKD_ID $ @;
            put MSKD_REF $ @;
            put MSKD_ALT $ @;
            put MSKD_QUAL @;
            put MSKD_FILTER $ @;
            put MSKD_INFO $ @;
            put MSKD_FORMAT $ @;
            put MSKD_GENOTYPE $ ;
            ;
        end;
        run;

%mend;
%iterdataset(dataset=design_file, function=%nrstr(%output_vcf(&line);));

/* Clean Up */
proc datasets nolist;
    delete out;
    delete line;
    delete w1118;
    run; quit;
