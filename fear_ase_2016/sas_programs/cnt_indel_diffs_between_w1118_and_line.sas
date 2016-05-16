/*******************************************************************************
* Filename: cnt_indel_diffs_between_w1118_and_line.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Merge w1118 with each line and count the number of snps they
* share. I think this step can be skipped, but need to check.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname THUMP '!HOME/thumper/cegs_ase_paper/indels';
libname dmel '!MCLAB/useful_dmel_data/flybase557/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/macros';
options SASAUTOS=(sasautos mymacros);
*/

/*
data design_file;
    set CEGS.genotype_list;
    run;
*/

%macro count_w11182line_snps(ID);

    proc sort data = THUMP.w1118_indel_vcf;
        by chrom pos;
        run;

    proc sort data = THUMP.&ID._indel_vcf;
        by chrom pos;
        run;

    data w1118_2_&ID._both  w1118_2_&ID._not_both  ;
        merge THUMP.w1118_indel_vcf (in=in1)  THUMP.&ID._indel_vcf (in=in2);
        by chrom pos;
        if in1 and in2 then output w1118_2_&ID._both ;
        else output w1118_2_&ID._not_both ;
        %*else if in1=1 and in2=0 then output &ID._only ;
        %*else if in1=0 and in2=1 then output w1118_&ID._only ;
        run;

    data w1118_2_&ID._both2 ;
        set w1118_2_&ID._both ;
        if Alt_&ID = ' ' then alt2_&ID = ref_&ID ;   /* if line snp is missing then is same as ref base */
        else alt2_&ID = alt_&ID ;
        if Alt_w1118 = ' ' then alt2_w1118 = ref_w1118 ;  /* if w1118 snp is missing then is same as ref base */
        else alt2_w1118 = alt_w1118 ;
        run ;

    data w1118_2_&ID._cnts;
        retain chrom pos snp_present w1118_snp &ID._snp ref_w1118 ref_&ID alt2_&ID 
		alt2_w1118 alt_&ID alt_w1118 ;
        set w1118_2_&ID._both2 ;

        if alt2_w1118 ne ref_w1118 then w1118_snp = 1 ;
        else w1118_snp = 0;
        if alt2_&ID ne ref_w1118 then &ID._snp = 1 ;
        else &ID._snp = 0 ;

        w1118_cnt = count(w1118, 'Het');
        if w1118_cnt > 0 then w1118_ishet = 1 ;
        else w1118_ishet = 0 ;
        &ID._cnt = count(&ID, 'Het');
        if &ID._cnt > 0 then &ID._ishet = 1 ;
        else &ID._ishet = 0;

        if w1118_ishet = 1 or &ID._ishet = 1 then snp_present = 0;  /* if w1118 is a het or line is a het then snp not present */
        else if alt2_w1118 ne alt2_&ID then snp_present = 1 ;       /* if w1118 snp is not the same as the line snp then snp is present */
        if w1118_snp = 1 then snp_present = 1 ;                     /* if there is a w1118 snp then the snp is present  */
        if &ID._snp = 1 then snp_present = 1 ;                      /* if there is a line snp then the snp is present */
        else snp_present = 0;
        drop w1118_cnt &ID._cnt;
        run ;

    proc summary data = w1118_2_&ID._cnts;
        var snp_present ;
        output out = sum_snp_diffs_&ID sum =sum_snps_w1118_to_line ;
        run ;

    data sum_snp_diffs_w1118_&ID ;
        set sum_snp_diffs_&ID ;
        rename _freq_ = num_in_w1118_and_line;
        line = "&ID" ;
        drop _type_  ;
        run;

%mend;

%iterdataset(dataset=design_file, function=%nrstr(%count_w11182line_snps(&line);));


data THUMP.indel_diffs_w1118_2_line_completes ;
    retain line sum_snps_w1118_to_line num_in_w1118_and_line;
    set sum_snp_diffs_w1118_: ;
    drop num_in_w1118_and_line ;
    run;

/*
proc export data = THUMP.indel_diffs_w1118_2_line
outfile='/home/ammorse/mclab/McIntyre_Lab/cegs_sergey/reports/ase/snp_diffs_w1118_2_line_
completes.csv'
dbms=csv replace ;
run;
*/


/* Clean Up */
proc datasets nolist;
    delete sum_snp_diffs_r101;
    delete sum_snp_diffs_r109;
    delete sum_snp_diffs_w1118_r101;
    delete sum_snp_diffs_w1118_r109;
    delete w1118_2_r101_both;
    delete w1118_2_r101_both2;
    delete w1118_2_r101_cnts;
    delete w1118_2_r101_not_both;
    delete w1118_2_r109_both;
    delete w1118_2_r109_both2;
    delete w1118_2_r109_cnts;
    delete w1118_2_r109_not_both;
    run; quit;

