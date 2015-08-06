/*******************************************************************************
* Filename: splicing_build_dsx_table_to_test_splicing.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: Build a table containing the InR fusions that are present in the
* cegsV data.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname cegs '!MCLAB/svn_lmm_dros_head_data/mel_cegs_expression_1st_106_F/OE_normalization/sas_data';
libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Create a list of InR Fusions that are in the cegsV data */
    data cegs_inr;
        set SEM.InR_gene_and_fusion;
        drop line FBgn0013984 Sxl dsx;
        run;

    proc transpose data=cegs_inr out=cegs_inr2;
        var :;
        run;

    data cegs_inr3;
        set cegs_inr2;
        label _name_ = ' ';
        rename _name_ = fusion_id;
        keep _name_;
        run;

/* Grab those fusions from the DSX data */
    proc sort data=cegs_inr3;
        by fusion_id;
        run;

    proc sort data=DSX.dsx_all;
        by fusion_id;
        run;

    data merged;
        retain fusion_id trt trt2 rep logrpkm;
        merge DSX.dsx_all (in=in1) cegs_inr3 (in=in2);
        by fusion_id;
        if in1 and in2;
        if trt eq 'AH_BerF' or trt eq 'AH_CSFemale' or trt eq 'AH_dsxNullF';
        if trt eq 'AH_BerF' or trt eq 'AH_CSFemale' then trt2 = 'control';
        else if trt eq 'AH_dsxNullF' then trt2='null';
        keep fusion_id trt trt2 rep logrpkm;
        run;

    proc sort data=merged;
        by fusion_id trt rep;
        run;

    data SEM.dsxnull_inr_fusions;
        set merged;
        run;

/* Clean Up */
proc datasets ;
    delete cegs_inr;
    delete cegs_inr2;
    delete cegs_inr3;
    delete merged;
    run; quit;
