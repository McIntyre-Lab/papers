/*******************************************************************************
* Filename: misc_test_cegsV_inr_raw_table.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Build table of the raw counts for the inr fusions. Add on the
* overall gene value that I am using for SEMs.
*******************************************************************************/

/* Libraries
    libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
    libname cegs '!MCLAB//svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
    libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
    filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
    options SASAUTOS=(sasautos mymacros);
*/

/* What Fusions are in InR */
    data inr_fus;
        set DMEL551.fb551_si_fusions;
        where exon_gene_id ? 'FBgn0013984';
        len = end - start;
        keep fusion_id len start;
        run;

    * The dmel table is unique by exon, so 'Fusions' are not unique;
    proc sort data=inr_fus nodupkey;
        by fusion_id;
        run;

/* Figure out which InR fusions are dropped during normalization */
    proc sort data=CEGS.flag_drop_fusion_by_line;
        by fusion_id;
        run;

    data flag_drop;
        merge CEGS.flag_drop_fusion_by_line (where=(mating_status eq 'V') in=in1) inr_fus (in=in2);
        by fusion_id;
        if in2;
        run;

    proc sort data=flag_drop;
        by start;
        run;

    proc export data=flag_drop outfile='!MCLAB/cegs_sem_sd_paper/reports/inr_fusions.csv' dbms=csv replace;
        putnames=yes;
        run;

/* Grab raw counts for InR fusions */
    * Lauren and I went through the InR wiggle plots and found that these three
    * fusions look to vary the most consistently between the different lines.
    ;
    proc sort data=CEGS.line_norm_centered;
        by fusion_id;
        run;

    data raw;
        merge CEGS.line_norm_centered (where=(mating_status eq 'V') in=in1) inr_fus (in=in2);
        by fusion_id;
        if in1 and in2;
        run;

    data raw2;
        retain line fusion_id apn uq_log_uq_center;
        set raw;
        keep line fusion_id apn uq_log_uq_center;
        run;

    proc sort data=raw2;
        by line;
        run;

    proc transpose data=raw2 out=apnflip;
        by line;
        var apn;
        id fusion_id;
        run;

    data apnflip2;
        set apnflip;
        drop _name_;
        run;

/* Grab normalized gene value for InR, dsx, and Sxl */
    data gene;
        set SEM.cegsV_by_gene_sbs;
        keep line FBgn0013984 dsx Sxl;
        run;

/* Combine Raw and gene level data */
    proc sort data=gene;
        by line;
        run;

    proc sort data=apnflip2;
        by line;
        run;

    data merged;
        merge gene (in=in1) apnflip2 (in=in2);
        by line;
        if in1;
        run;

    data SEM.inr_gene_and_fusion;
        set merged;
        run;

/* Clean up */
proc datasets ;
    delete apnflip;
    delete apnflip2;
    delete flag_drop;
    delete gene;
    delete inr_fus;
    delete merged;
    delete raw;
    delete raw2;
    run; quit;

