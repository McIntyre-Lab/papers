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

%macro output_vcf(ID) ;
    data w1118;
        retain chrom pos ID_w1118 REF_w1118 ALT_w1118 QUAL_w1118 FILTER_w1118 INFO_w1118 FORMAT_w1118 w1118;
        set THUMP.w1118_2_&ID._vcf_flags2;
        where alt_w1118 ne '' and Output_Lvl1_MSKD = 1;
        label chrom = #CHROM;
        label pos = POS;
        label ID_w1118 = ID;
        label REF_w1118 = REF;
        label ALT_w1118 = ALT;
        label QUAL_w1118 = QUAL;
        label FILTER_w1118 = FILTER;
        label INFO_w1118 = INFO;
        label FORMAT_w1118 = FORMAT;
        label w1118 = GENOTYPE;
        keep chrom pos ID_w1118 REF_w1118 ALT_w1118 QUAL_w1118 FILTER_w1118 INFO_w1118 FORMAT_w1118 w1118;
        run;

	data _null_;
		file "!HOME/sandbox/cegs_ase_paper/ase_lvl1_filtered_vcf_files/w1118_w11182&ID._lvl1.vcf"
		delimiter='09'x DSD DROPOVER lrecl=32767;
		if _n_ = 1 then do;
			put '#CHROM' '09'x 'POS' '09'x 'ID' '09'x 'REF' '09'x 'ALT' '09'x 'QUAL' '09'x 
				'FILTER' '09'x 'INFO' '09'x 'FORMAT' '09'x 'GENOTYPE';
		end;
		set  w1118   end=EFIEOD;
		format chrom $32. ;
		format pos best12. ;
		format ID_w1118 best12. ;
		format REF_w1118 $300. ;
		format ALT_w1118 $300. ;
		format QUAL_w1118 best12. ;
		format FILTER_w1118 $36. ;
		format INFO_w1118 $262. ;
		format FORMAT_w1118 $18. ;
		format w1118 $74. ;
		do;
			put chrom $ @;
			put pos @;
			put ID_w1118 $ @;
			put REF_w1118$ @;
			put ALT_w1118 $ @;
			put QUAL_w1118 @;
			put FILTER_w1118 $ @;
			put INFO_w1118 $ @;
			put FORMAT_w1118 $ @;
			put w1118 $ ;
			;
		end;
		run;

    data line;
        retain chrom pos ID_&ID. REF_&ID. ALT_&ID. QUAL_&ID. FILTER_&ID. INFO_&ID. FORMAT_&ID. &ID.;
        set THUMP.w1118_2_&ID._vcf_flags2;
        where alt_&ID. ne '' and Output_Lvl1_MSKD = 1;
        label chrom = #CHROM;
        label pos = POS;
        label ID_&ID. = ID;
        label REF_&ID. = REF;
        label ALT_&ID. = ALT;
        label QUAL_&ID. = QUAL;
        label FILTER_&ID. = FILTER;
        label INFO_&ID. = INFO;
        label FORMAT_&ID. = FORMAT;
        label &ID. = GENOTYPE;
        keep chrom pos ID_&ID. REF_&ID. ALT_&ID. QUAL_&ID. FILTER_&ID. INFO_&ID. FORMAT_&ID. &ID.;
        run;

	data _null_;
		file "!HOME/sandbox/cegs_ase_paper/ase_lvl1_filtered_vcf_files/&ID._w11182&ID._lvl1.vcf"
		delimiter='09'x DSD DROPOVER lrecl=32767;
		if _n_ = 1 then do;
			put '#CHROM' '09'x 'POS' '09'x 'ID' '09'x 'REF' '09'x 'ALT' '09'x 'QUAL' '09'x 
				'FILTER' '09'x 'INFO' '09'x 'FORMAT' '09'x 'GENOTYPE';
		end;
		set  line   end=EFIEOD;
		format chrom $32. ;
		format pos best12. ;
		format ID_&ID. best12. ;
		format REF_&ID. $300. ;
		format ALT_&ID. $300. ;
		format QUAL_&ID. best12. ;
		format FILTER_&ID. $36. ;
		format INFO_&ID. $262. ;
		format FORMAT_&ID. $18. ;
		format &ID. $74. ;
		do;
			put chrom $ @;
			put pos @;
			put ID_&ID. $ @;
			put REF_&ID.$ @;
			put ALT_&ID. $ @;
			put QUAL_&ID. @;
			put FILTER_&ID. $ @;
			put INFO_&ID. $ @;
			put FORMAT_&ID. $ @;
			put &ID. $ ;
			;
		end;
		run;
%mend;
%iterdataset(dataset=design_file, function=%nrstr(%output_vcf(&line);));

proc datasets nolist;
    delete w1118;
    delete line;
    run; quit;
