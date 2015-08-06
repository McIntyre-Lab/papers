/*******************************************************************************
* Filename: misc_test_build_table_to_test_splicing.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Create a stacked dataset with sex det and InR genes.
*
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname cegs '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* Pull virgin out of normalized and centered data */
    data virgin_w_rep;
        set CEGS.ccfus_norm_centered;
        where mating_status = 'V';
        keep fusion_id line rep uq_log_uq_center;
        run;

    proc sort data=virgin_w_rep;
        by fusion_id line;
        run;

    proc sort data=SEM.cegs_virgin_norm_cent;
        by fusion_id line;
        run;

    data merged;
        merge virgin_w_rep (in=in1) SEM.cegs_virgin_norm_cent (in=in2);
        by fusion_id line;
        if in2;
        keep fusion_id line rep mean_exp uq_log_uq_center symbol_cat;
        run;


/* Grab the Sex Det genes */
    proc sort data=merged;
        by fusion_id;
        run;

    proc sort data=SEM.cegs_flag_sex_det;
        by fusion_id;
        run;

    data norm_sex;
        merge merged (in=in1) SEM.cegs_flag_sex_det (in=in2);
        by fusion_id;
        if flag_sex_det = 1;
        drop flag_sex_det;
        run;

    * clean up cat'd gene symbols;
        proc freq data=norm_sex;
            table symbol_cat;
            run;

        data norm_sex2;
            set norm_sex;
            if symbol_cat eq 'CG10841|sqd' then symbol_cat = 'sqd';
            if symbol_cat eq 'CG34424|vir' then symbol_cat = 'vir';
            if symbol_cat eq 'Dek|Psi' then symbol_cat = 'Psi';
            if symbol_cat eq 'blos1|tra2' then symbol_cat = 'tra2';
            if symbol_cat eq 'dsx|lds' then symbol_cat = 'dsx';
            if symbol_cat eq 'fl(2)d' then symbol_cat = 'fl_2_d';
            if symbol_cat eq 'l(3)73Ah|tra' then symbol_cat = 'tra';
            run;

/* Grab the gene in the InR region */
    data _null_;
        set DMEL551.symbol2coord;
        where symbol eq 'InR';
        call symput('chrom', chrom);
        call symput('start', start);
        call symput('end', end);
        run;

    data fusList;
        set DMEL551.fb551_si_fusions_unique_flagged;
        where chrom eq "&chrom" and start ge &start and end le &end;
        keep fusion_id;
        run; *18 obs including CR43653;

    proc sort data=fusList;
        by fusion_id;
        run;

    data norm_inr;
        merge merged (in=in1) fusList (in=in2);
        by fusion_id;
        if in2;
        if uq_log_uq_center = '.' then delete; * Delete fusions that don't have any expression;
        run;

/* Combine datasets and make permanent splice dataset */
    data SEM.cegsV_splice_data;
        set norm_sex2 norm_inr;
        run;

/* Clean Up */
proc datasets nolist;
    delete virgin_w_rep;
    delete merged;
    delete fuslist;
    delete norm_inr;
    delete norm_sex;
    delete norm_sex2;
    run; quit;
