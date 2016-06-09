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
libname THUMP '!HOME/thumper/cegs_ase_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase557/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/********************************************************************************/
/* Categorize SNP classes                                                       */
/* Flag SNPS for Level1 outputting                                              */
/* Output VCF files for creating masked alignments                              */
/* do not delete THUMP.w1118_2_&ID._vcf_flags2 !!                               */
/********************************************************************************/

/*
data design_file;
    set CEGS.genotype_list;
    run;
*/

%macro flag_w11182line_vcf (ID);
    data THUMP.w1118_2_&ID._vcf_flags2;
        set WORK.w1118_2_&ID._vcf;
        length MSKD_REF $ 300;
        length MSKD_ALT $ 300;
        length MSKD_FILTER $ 36;
        length MSKD_INFO $ 262;
        length MSKD_FORMAT $ 18;
        length MSKD_ID  $ 26; ;
        length MSKD_GENOTYPE $ 74;

        *Make sure that if the SNP is called in both, the reference base given is the same in both vcf files shouldnt be any ref_not_identical = 1;
        if VCF_in_both = 1 and REF_w1118 ~= REF_&ID then REF_NOT_IDENTICAL = 1;
        else REF_NOT_IDENTICAL = 0;
        if VCF_in_both ~= 1 then REF_NOT_IDENTICAL = ".";

        *Case 1, either position is called as het;
        if w1118_ishet = 1 or &ID._ishet = 1 then do;
            Output_Lvl1_MSKD = 0;
            Output_Lvl1_UPD = 0;
            case = 1;
        end;

        *Case 2, both are the reference;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 = ALT_&ID and ALT_w1118 = "" and ALT_&ID = "" then do;
            Output_Lvl1_MSKD = 0;
            Output_Lvl1_UPD = 0;
            case = 2;
        end;

        *Case 3, both have a SNP as compared to FB5.30 and the base at the SNP position is the same in both w1118 and the line;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 = ALT_&ID and ALT_w1118 ~= "" and ALT_&ID ~= "" then do;
            MSKD_ID = ID_&ID;
            MSKD_REF = REF_&ID;
            MSKD_ALT = ALT_&ID; 
            MSKD_QUAL = QUAL_&ID;
            MSKD_FILTER = FILTER_&ID ;
            MSKD_INFO = INFO_&ID;
            MSKD_FORMAT = FORMAT_&ID;
            MSKD_GENOTYPE = &ID.;

            *If you want to filter these positions before masking replace the current flags with the below.
            *AND remove the above lines that assign MSKD vcf variables;
            *I defaulted to outputing the line here, but if they are the same it should not matter.
            *Only the position is used in masked updating;

            *Output_Lvl1_MSKD = 0;
            *Output_Lvl1_UPD = 1;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 3;
        end;

        *Case 4, only ALT_w1118 or ALT_&ID is the reference;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 ~= ALT_&ID and ALT_w1118 = "" then do;
            MSKD_ID = ID_&ID;
            MSKD_REF = REF_&ID;
            MSKD_ALT = ALT_&ID; 
            MSKD_QUAL = QUAL_&ID;
            MSKD_FILTER = FILTER_&ID ;
            MSKD_INFO = INFO_&ID;
            MSKD_FORMAT = FORMAT_&ID;
            MSKD_GENOTYPE = &ID.;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 4;
        end;

        if w1118_ishet = 0 and VCF_in_both = 1 and ALT_w1118 ~= ALT_&ID and ALT_&ID = "" then do;
            MSKD_ID = ID_w1118;
            MSKD_REF = REF_w1118; 
            MSKD_ALT = ALT_w1118 ;
            MSKD_QUAL = QUAL_w1118;
            MSKD_FILTER = FILTER_w1118; 
            MSKD_INFO = INFO_w1118;
            MSKD_FORMAT = FORMAT_w1118; 
            MSKD_GENOTYPE = w1118 ;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 4;
        end;

        *Case 5, where both are a SNP from FB5.30, but each is a different base;
        *The information output to the VCF file for masking defaults to the base for the line;
        *However, only the position information is used at this point;
        if w1118_ishet = 0 and &ID._ishet = 0 and VCF_in_both = 1 and ALT_w1118 ~= ALT_&ID and ALT_&ID ~= "" and ALT_w1118 ~= "" then do;
            MSKD_ID = ID_&ID;
            MSKD_REF = REF_&ID; 
            MSKD_ALT = ALT_&ID ;
            MSKD_QUAL = QUAL_&ID;
            MSKD_FILTER = FILTER_&ID; 
            MSKD_INFO = INFO_&ID;
            MSKD_FORMAT = FORMAT_&ID; 
            MSKD_GENOTYPE = &ID.;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 5;
        end;

        *Case 6, the SNP is found in the ID vcf file only or in the w1118 vcf file only and is not a het;
        if &ID._ishet = 0 and VCF_in_&ID._only = 1 and ALT_&ID ~= "" then do;
            MSKD_ID = ID_&ID;
            MSKD_REF = REF_&ID; 
            MSKD_ALT = ALT_&ID; 
            MSKD_QUAL = QUAL_&ID;
            MSKD_FILTER = FILTER_&ID; 
            MSKD_INFO = INFO_&ID;
            MSKD_FORMAT = FORMAT_&ID;
            MSKD_GENOTYPE = &ID. ;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 6;
        end;

        if w1118_ishet = 0 and VCF_in_w1118_only = 1 and ALT_w1118 ~= "" then do;
            MSKD_ID = ID_w1118;
            MSKD_REF = REF_w1118;
            MSKD_ALT = ALT_w1118; 
            MSKD_QUAL = QUAL_w1118;
            MSKD_FILTER = FILTER_w1118; 
            MSKD_INFO = INFO_w1118;
            MSKD_FORMAT = FORMAT_w1118;
            MSKD_GENOTYPE = w1118 ;
            Output_Lvl1_MSKD = 1;
            Output_Lvl1_UPD = 1;
            case = 6;
        end;

        *Case 7, the position is found in the ID vcf file only or in the w1118 vcf file only and is not a het, but it is called as the reference;
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
       delete w1118_2_&ID._vcf;
       run;
%mend;
%iterdataset(dataset=design_file, function=%nrstr(%flag_w11182line_vcf(&line);));

 *Output a VCF file for masking ;
%macro output_vcf(ID) ;
    data w11182&ID._vcf_MSKD_out1;
        set THUMP.w1118_2_&ID._vcf_flags2;
        if Output_Lvl1_MSKD = 1;
        label chrom = #CHROM;
        label pos = POS;
        label MSKD_ID = ID;
        label MSKD_REF = REF;
        label MSKD_ALT = ALT;
        label MSKD_QUAL = QUAL;
        label MSKD_FILTER = FILTER;
        label MSKD_INFO = INFO;
        label MSKD_FORMAT = FORMAT;
        label MSKD_GENOTYPE = GENOTYPE;
        keep chrom pos MSKD_ID MSKD_REF MSKD_ALT MSKD_QUAL MSKD_FILTER MSKD_INFO MSKD_FORMAT MSKD_GENOTYPE;
        run;

    data w11182&ID._vcf_MSKD_out;
        retain chrom pos MSKD_ID MSKD_REF MSKD_ALT MSKD_QUAL MSKD_FILTER MSKD_INFO MSKD_FORMAT MSKD_GENOTYPE;
        set w11182&ID._vcf_MSKD_out1;
        run;

	data _null_;
		file "!HOME/sandbox/cegs_ase_paper/ase_lvl1_filtered_vcf_files/w11182&ID._MSKD.vcf"
		delimiter='09'x DSD DROPOVER lrecl=32767;
		if _n_ = 1 then do;
			put '#CHROM' '09'x 'POS' '09'x 'ID' '09'x 'REF' '09'x 'ALT' '09'x 'QUAL' '09'x 
				'FILTER' '09'x 'INFO' '09'x 'FORMAT' '09'x 'GENOTYPE';
		end;
		set  W11182&ID._VCF_MSKD_OUT   end=EFIEOD;
		format chrom $32. ;
		format pos best12. ;
		format MSKD_ID $26. ;
		format MSKD_REF $300. ;
		format MSKD_ALT $300. ;
		format MSKD_QUAL best12. ;
		format MSKD_FILTER $36. ;
		format MSKD_INFO $262. ;
		format MSKD_FORMAT $18. ;
		format MSKD_GENOTYPE $74. ;
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

    proc datasets nolist;
        delete W11182&ID._VCF_MSKD_OUT;
        delete W11182&ID._VCF_MSKD_OUT1;
        run; quit;

%mend;
%iterdataset(dataset=design_file, function=%nrstr(%output_vcf(&line);));
