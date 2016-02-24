/*******************************************************************************
* Filename: level2_export_UPD_vcf.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Split VCFs for w1118 and line, keeping the SNPs and indels
* together. Then export a vcf for updating.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

%macro split_UPD(ID);
    /* Determine if the alt is a SNP or an indel and make output for Line */
    data &ID._out;
        set THUMP.flag_lvl2_w1118_2_&ID;
        where line eq "&ID" and OUTPUT_Lvl2_UPD eq 1;
        keep chrom pos UPD_ID UPD_REF UPD_ALT UPD_QUAL UPD_FILTER UPD_INFO UPD_FORMAT UPD_GENOTYPE;
        run;

    /* Determine if the alt is a SNP or an indel and make output for w1118 */
    data w1118_out;
        set THUMP.flag_lvl2_w1118_2_&ID;
        where line eq "w1118" and OUTPUT_Lvl2_UPD eq 1;
        keep chrom pos UPD_ID UPD_REF UPD_ALT UPD_QUAL UPD_FILTER UPD_INFO UPD_FORMAT UPD_GENOTYPE;
        run;

    /* Create VCF Files */
    data _null_;
        file "!HOME/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files/&ID._w11182&ID._UPD.vcf"
        delimiter='09'x DSD DROPOVER lrecl=32767;
        if _n_ = 1 then do;
            put '#CHROM' '09'x 'POS' '09'x 'ID' '09'x 'REF' '09'x 'ALT' '09'x 'QUAL' '09'x 
                'FILTER' '09'x 'INFO' '09'x 'FORMAT' '09'x 'GENOTYPE';
        end;
        set  &ID._out   end=EFIEOD;
        format chrom $32. ;
        format pos best12. ;
        format UPD_ID best12. ;
        format UPD_REF $300. ;
        format UPD_ALT $300. ;
        format UPD_QUAL best32. ;
        format UPD_FILTER $36. ;
        format UPD_INFO $262. ;
        format UPD_FORMAT $18. ;
        format UPD_GENOTYPE $300. ;
        do;
            put chrom $ @;
            put pos @;
            put UPD_ID $ @;
            put UPD_REF $ @;
            put UPD_ALT $ @;
            put UPD_QUAL @;
            put UPD_FILTER $ @;
            put UPD_INFO $ @;
            put UPD_FORMAT $ @;
            put UPD_GENOTYPE $ ;
            ;
        end;
        run;

    data _null_;
        file "!HOME/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files/w1118_w11182&ID._UPD.vcf"
        delimiter='09'x DSD DROPOVER lrecl=32767;
        if _n_ = 1 then do;
            put '#CHROM' '09'x 'POS' '09'x 'ID' '09'x 'REF' '09'x 'ALT' '09'x 'QUAL' '09'x 
                'FILTER' '09'x 'INFO' '09'x 'FORMAT' '09'x 'GENOTYPE';
        end;
        set  w1118_out   end=EFIEOD;
        format chrom $32. ;
        format pos best12. ;
        format UPD_ID best12. ;
        format UPD_REF $300. ;
        format UPD_ALT $300. ;
        format UPD_QUAL best32. ;
        format UPD_FILTER $36. ;
        format UPD_INFO $262. ;
        format UPD_FORMAT $18. ;
        format UPD_GENOTYPE $300. ;
        do;
            put chrom $ @;
            put pos @;
            put UPD_ID $ @;
            put UPD_REF $ @;
            put UPD_ALT $ @;
            put UPD_QUAL @;
            put UPD_FILTER $ @;
            put UPD_INFO $ @;
            put UPD_FORMAT $ @;
            put UPD_GENOTYPE $ ;
            ;
        end;
        run;

    /* Export List of regions to Mask in BED format */
        data _null_;
            file "!HOME/sandbox/cegs_ase_paper/ase_lvl2_filtered_vcf_files/pos_to_permMask_w11182&ID..bed"
            delimiter='09'x DSD DROPOVER lrecl=32767;
            set THUMP.pos_to_permmask_w11182&ID  end=EFIEOD;
            format chrom $32. ;
            format start best12. ;
            format end best12. ;
            do;
                put chrom $ @;
                put start @;
                put end $;
                ;
            end;
            run;

    /* Clean Up */
        proc datasets nolist;
            delete &ID._OUT;
            delete W1118_OUT;
            run; quit;
%mend;
%iterdataset(dataset=design_file, function=%nrstr(%split_UPD(&line);));

