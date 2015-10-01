/*******************************************************************************
* Filename: level2_merge_level1_filters.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Merge RNA counts to level 1 filters and make level 2 filter
* Lvl2 Filter = RNA CVG > 0 or DNA CVG > 5
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/*
data design_file;
    set CEGS.genotype_list;
    run;

    %let ID=w79;
*/

%macro merge_lvl1(ID);
    /* Stack Lvl1 filters for SNPS and Indels */
        data w1118_snp;
            length w1118 $ 300.;
            length line $ 5.;
            length type $ 5.;
            set THUMP.w1118_2_&ID._vcf_flags2;
            where output_lvl1_MSKD = 1 and ALT_w1118 ne '';
            rename REF_w1118        = UPD_REF;
            rename ALT_w1118        = UPD_ALT;
            rename FILTER_w1118     = UPD_FILTER;
            rename INFO_w1118       = UPD_INFO;
            rename QUAL_w1118       = UPD_QUAL;
            rename FORMAT_w1118     = UPD_FORMAT;
            rename w1118            = UPD_GENOTYPE;
            rename ID_w1118         = UPD_ID;
            line = 'w1118';
            type = 'snp';
            keep CHROM POS VCF_in_both VCF_in_w1118_only VCF_in_&ID._only w1118_ishet
            &ID._ishet REF_NOT_IDENTICAL Output_Lvl1_MSKD Output_Lvl1_UPD case
            REF_w1118 ALT_w1118 FILTER_w1118 INFO_w1118 QUAL_w1118 FORMAT_w1118
            w1118 ID_w1118 line type;
            run;

        data &ID._snp;
            length &ID. $ 300.;
            length line $ 5.;
            length type $ 5.;
            set THUMP.w1118_2_&ID._vcf_flags2;
            where output_lvl1_MSKD = 1 and ALT_&ID ne '';
            rename REF_&ID.        = UPD_REF;
            rename ALT_&ID.        = UPD_ALT;
            rename FILTER_&ID.     = UPD_FILTER;
            rename INFO_&ID.       = UPD_INFO;
            rename QUAL_&ID.       = UPD_QUAL;
            rename FORMAT_&ID.     = UPD_FORMAT;
            rename &ID.            = UPD_GENOTYPE;
            rename ID_&ID.         = UPD_ID;
            line = "&ID.";
            type = "snp";
            keep CHROM POS VCF_in_both VCF_in_w1118_only VCF_in_&ID._only w1118_ishet
            &ID._ishet REF_NOT_IDENTICAL Output_Lvl1_MSKD Output_Lvl1_UPD case
            REF_&ID. ALT_&ID. FILTER_&ID. INFO_&ID. QUAL_&ID. FORMAT_&ID.
            &ID. ID_&ID. line type;
            run;

        data w1118_indel;
            length line $ 5.;
            length type $ 5.;
            set THUMP.w1118_2_&ID._indel_vcf_flags2;
            where output_lvl1_MSKD = 1 and ALT_w1118 ne '';
            rename REF_w1118        = UPD_REF;
            rename ALT_w1118        = UPD_ALT;
            rename FILTER_w1118     = UPD_FILTER;
            rename INFO_w1118       = UPD_INFO;
            rename QUAL_w1118       = UPD_QUAL;
            rename FORMAT_w1118     = UPD_FORMAT;
            rename w1118            = UPD_GENOTYPE;
            rename ID_w1118         = UPD_ID;
            line = 'w1118';
            type = 'indel';
            keep CHROM POS VCF_in_both VCF_in_w1118_only VCF_in_&ID._only w1118_ishet
            &ID._ishet REF_NOT_IDENTICAL Output_Lvl1_MSKD Output_Lvl1_UPD case
            REF_w1118 ALT_w1118 FILTER_w1118 INFO_w1118 QUAL_w1118 FORMAT_w1118
            w1118 ID_w1118 line type;
            run;

        data &ID._indel;
            length line $ 5.;
            length type $ 5.;
            set THUMP.w1118_2_&ID._indel_vcf_flags2;
            where output_lvl1_MSKD = 1 and ALT_&ID ne '';
            rename REF_&ID.        = UPD_REF;
            rename ALT_&ID.        = UPD_ALT;
            rename FILTER_&ID.     = UPD_FILTER;
            rename INFO_&ID.       = UPD_INFO;
            rename QUAL_&ID.       = UPD_QUAL;
            rename FORMAT_&ID.     = UPD_FORMAT;
            rename &ID.            = UPD_GENOTYPE;
            rename ID_&ID.         = UPD_ID;
            line = "&ID.";
            type = "indel";
            keep CHROM POS VCF_in_both VCF_in_w1118_only VCF_in_&ID._only w1118_ishet
            &ID._ishet REF_NOT_IDENTICAL Output_Lvl1_MSKD Output_Lvl1_UPD case
            REF_&ID. ALT_&ID. FILTER_&ID. INFO_&ID. QUAL_&ID. FORMAT_&ID.
            &ID. ID_&ID. line type;
            run;

        data variants;
            set w1118_snp w1118_indel &ID._snp &ID._indel;
            lref = length(UPD_REF);
            lalt = length(UPD_ALT);
            end = POS + lref - 1;
            run;

    /* Look for overlapping SNPs and InDels and flag them */
        proc sort data=variants;
            by chrom POS end;
            run;

        data find_overlap;
            set variants;
            by chrom pos end;
            informat curref $300.;
            informat curalt $300.;
            informat preref $300.;
            informat prealt $300.;
            informat preline $5.;
            informat curline $5.;
            /* Initilize values */
            retain regionstart regionend 
                   minMask maxMask
                   count 
                   flag_merge
                   curref curalt
                   preref prealt
                   preline curline
                   ;
            if first.chrom then do;
                regionend = .;
                regionstart = .;
                minMask = .;
                maxMask = .;
                count = .;
                flag_merge = .;
                curref = '';
                curalt = '';
                preref = '';
                prealt = '';
                preline = '';
                curline = '';
            end;
            /* Start a new set if pos is greater than the current end */
            if pos > regionend then do;
                if not first.chrom then do;
                    if count > 2 then flag_merge = 1;
                    else if count = 2 then do;
                        /* If there are two overlapping variants, but they are
                        * in different lines, then do not merge */
                        if (curref ne preref) or (curalt ne prealt) then do;
                            if preline ne curline then flag_merge = 0;
                            else flag_merge = 1;
                        end;
                    end;
                    output;
                end;
                regionstart=pos;
                regionend=end;
                maxMask = regionstart + lref - 1;
                minMask = regionstart;
                count = 1;
                flag_merge = 0;
                preref = curref;
                prealt = curalt;
                preline = curline;
                curline = line;
                curref = UPD_REF;
                curalt = UPD_ALT;
            end;
            /* Extend the set if pos is between regionstart and regionend, but pos + lref is past region end */
            else if end > regionend then do;
                regionend=end;
                if pos + lref - 1 > maxMask then maxMask = pos + lref - 1;
                count = count + 1;
                preref = curref;
                prealt = curalt;
                preline = curline;
                curline = line;
                curref = UPD_REF;
                curalt = UPD_ALT;
            end;
            else do;
                count = count + 1;
                preref = curref;
                prealt = curalt;
                preline = curline;
                curline = line;
                curref = UPD_REF;
                curalt = UPD_ALT;
            end;
            /* If last item in the chrom then output */
            if last.chrom then do;
                if count > 2 then flag_merge = 1;
                else if count = 2 then do;
                    if (curref ne preref) or (curalt ne prealt) then do;
                        if preline ne curline then flag_merge = 0;
                        else flag_merge = 1;
                    end;
                end;
                output;
            end;
            keep chrom maxMask minMask flag_merge;
            run;

        data overlapping;
            set find_overlap;
            where flag_merge eq 1;
            run;

        /* create permanent dataset */
            data THUMP.pos_to_permmask_w11182&ID;
                set overlapping;
                rename minmask = start;
                rename maxmask = end;
                keep chrom minmask maxmask;
                run;

        /* Create an overlap flag by expanding position ranges */
        %macro generate_overlap(chrom, start, end);
            data coords;
                informat chrom $32.;
                informat pos best32.;
                retain chrom pos;
                chrom = "&chrom";
                do pos = &start to &end;
                    output;
                end;
                run;

            proc append base=overlap_pos data=coords;
            run;
        %mend;
        %iterdataset(dataset=overlapping, function=%nrstr(%generate_overlap(&chrom, &minmask, &maxmask);));

        proc sort data=overlap_pos;
            by chrom pos;
            run;

    /* Merge on overlap flag to stack of variants*/
        proc sort data=variants;
            by chrom pos;
            run;

        data mask_perm;
            merge variants (in=in1) overlap_pos (in=in2);
            by chrom pos;
            if in2 then flag_overlap = 1;
            else flag_overlap = 0;
            if in1 then output;
            run;

    /* Merge RNA support counts */
        proc sort data=mask_perm ;
            by chrom pos UPD_ALT;
            run;

        proc sort data=&ID._counts;
            by chrom pos alt; 
            run;

        data merged;
            merge mask_perm (in=in1) &ID._counts (in=in2 rename=(alt=UPD_ALT));
            by chrom pos UPD_ALT;
            if in1;
            run; /* This merge gives a duplicate key warning because line is not included in the by statement, it is ok */

        proc sort data=merged ;
            by chrom pos;
            run;

        /* Flag SNPs and Indels that have RNA support in both line and tester as verified */
            proc transpose data=merged out=both delimiter=_;
                by chrom pos;
                var altRNACnt;
                id line type;
                run;

            data flag_verified;
                set both;
                if &ID._snp > 0 and w1118_snp > 0 then both_snp_verified = 1;
                else both_snp_verified = 0;
                if &ID._indel > 0 and w1118_indel > 0 then both_indel_verified = 1;
                else both_indel_verified = 0;
                keep chrom pos both_snp_verified both_indel_verified;
                run;

                /* Quick check
                    proc freq data=flag_verified;
                        table both_snp_verified; * 24554 with both;
                        table both_indel_verified; * 190 with both;
                        run;
                */
                
        /* Create Level 2 Flag, where (RNA support must be greater than 0 OR
         * DNA support is greater than 5) and variant does not overlap with
         * other variants */
        data THUMP.flag_lvl2_w1118_2_&ID.;
            merge merged (in=in1) flag_verified (in=in2);
            by chrom pos ;
            if altRnaCnt ge 10 then flag_RNA_10X_cvg = 1;
            else flag_RNA_10X_cvg  = 0;
            if DnaCnt ge 5 then flag_DNA_cvg_5 = 1;
            else flag_DNA_cvg_5  = 0;
            if (altRNACnt > 0 or flag_DNA_cvg_5 eq 1) and flag_overlap = 0 then Output_Lvl2_UPD = 1;
            else Output_Lvl2_UPD = 0;
            run;

    /* Clean up */
    proc datasets nolist;
        delete both;
        delete flag_verified;
        delete merged;
        delete &ID._indel;
        delete &ID._snp;
        delete &ID._counts;
        delete variants;
        delete w1118_indel;
        delete w1118_snp;
        delete find_overlap;
        delete mask_perm;
        delete overlapping;
        delete overlap_pos;
        delete coords;
        run; quit;
%mend;
%iterdataset(dataset=design_file, function=%nrstr(%merge_lvl1(&line);));
