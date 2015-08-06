/*******************************************************************************
* Filename: program_project_upstream_fru_dsx_not_ds_sxl.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dmel '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Grab models of interest and flag if downstream of sxl */
    data ag;
        set SEM.cegsv_ag_w_bic_cutoff;
        where model eq 'Model_15' or model eq 'Model_16' or model eq 'Model_26' or
        model eq 'Model_27' or model eq 'Model_34' or model eq 'Model_35' or model
        eq 'Model_1' ;
        if flag_expanded eq 1;
        if model eq 'Model_1' then flag_ds_sxl = 1;
        else flag_ds_sxl = 0;
        keep primary_fbgn symbol model flag_ds_sxl;
        run;

/* Only keep models that are upstream of dsx and fru and not downstream of sxl */
    proc sort data=ag;
        by primary_fbgn;
        run;

    proc freq data=ag noprint;
        by primary_fbgn;
        tables flag_ds_sxl /out=freqs;
        run;

    data us_dsx_fru;
        set freqs;
        if flag_ds_sxl eq 0 and percent eq 100;
        keep primary_fbgn;
        run;

/* Merge on gene symbol */
    proc sort data=dmel.symbol2fbgn;
        by primary_fbgn;
        run;

    data merged;
        merge us_dsx_fru (in=in1) dmel.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1 then flag_us_fru_dsx = 1;
        else flag_us_fru_dsx = 0;
        run;

/* Merge to ribo 575 gene list */
    libname ribo '!MCLAB/arbeitman/arbeitman_ribotag/sas_data';

    * 5 FBGNs did not merge cleanely, renaming by hand;
    data ribo;
        set RIBO.intersectipanddalton_575;
        if primary_fbgn eq 'FBgn0005634' then primary_fbgn = 'FBgn0265434';
        if primary_fbgn eq 'FBgn0008654' then primary_fbgn = 'FBgn0265623';
        if primary_fbgn eq 'FBgn0261642' then primary_fbgn = 'FBgn0265487';
        if primary_fbgn eq 'FBgn0262110' then primary_fbgn = 'FBgn0265988';
        if primary_fbgn eq 'FBgn0263218' then primary_fbgn = 'FBgn0265296';
        keep primary_fbgn;
        run;

    proc sort data=ribo;
        by primary_fbgn;
        run;

    proc sort data=merged;
        by primary_fbgn;
        run;

    data ribomerge tocheck;
        merge  merged (in=in1) ribo (in=in2);
        by primary_fbgn;
        if in2 then flag_575 = 1;
        else flag_575 = 0;
        if in2 and not in1 then output tocheck;
        else output ribomerge;
        run;

    proc freq data=ribomerge;
        tables flag_us_fru_dsx*flag_575 /chisq exact expect;
        run;

    /*
        flag_us_fru_dsx     flag_575

        Frequency|
        Expected |       0|       1|  Total
        ---------+--------+--------+
               0 | 228631 |    571 | 229202
                 | 228627 | 574.71 |
        ---------+--------+--------+
               1 |    111 |      4 |    115
                 | 114.71 | 0.2884 |
        ---------+--------+--------+
        Total      228742      575   229317

        Statistic      DF       Value      Prob
        ---------------------------------------
        Chi-Square      1     47.9194    <.0001
    */

/* Merge on Fru BS and over expression */
    libname FRU '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';

    data fru;
        set FRU.motif_flags_and_cnts;
        keep primary_fbgn flag_fru_a_motif flag_fru_b_motif flag_fru_c_motif;
        run;

    data fru2;
        set FRU.flag_ind_rep;
        sums = sum(flag_male_ind, flag_male_rep, flag_female_ind, flag_female_rep);
        if sums > 0 then flag_fru_over = 1;
        else flag_fru_over = 0;
        keep primary_fbgn flag_fru_over;
        run;

    proc sort data=fru;
        by primary_fbgn;
        run;

    proc sort data=fru2;
        by primary_fbgn;
        run;

    data fru3;
        merge fru (in=in1) fru2 (in=in2);
        by primary_fbgn;
        if in1 and not in2 then flag_fru_over = 0;
        if in2 and not in1 then do;
            flag_fru_a_motif = 0;
            flag_fru_b_motif = 0;
            flag_fru_c_motif = 0;
        end;
        run;

    /* Merge on DSX Null expression */
        data dsx_null;
            set SEM.dsxnullf_induced SEM.dsxnullf_repressed;
            rename fbgn_cat = primary_fbgn;
            run;

        proc sort data=dsx_null nodupkey;
            by primary_fbgn;
            run;

        data dsxfru;
            merge fru3 (in=in1) dsx_null (in=in2);
            by primary_fbgn;
            if in2 then flag_dsx_null = 1;
            else flag_dsx_null = 0;
            run;

    /* Convert Fru FBgn from 5.30 to 5.51 */
        proc sort data=dsxfru;
            by primary_fbgn;
            run;

        proc sort data=dmel.fbgn2oldfbgn;
            by other_fbgn;
            run;

        data frufb551;
            retain primary_fbgn;
            merge dsxfru (in=in1 rename=(primary_fbgn = fb530_fbgn)) dmel.fbgn2oldfbgn (in=in2 rename=(other_fbgn = fb530_fbgn));
            by fb530_fbgn;
            if in1 and in2;
            run;

        proc sort data=frufb551;
            by primary_fbgn;
            run;

        proc means data=frufb551 noprint;
            by primary_fbgn;
            output out=sums sum(flag_fru_a_motif)= sum(flag_fru_b_motif)=
            sum(flag_fru_c_motif)= sum(flag_fru_over)= sum(flag_dsx_null)=
            /autoname;
            run;

        data frufb551;
            set sums;
            if flag_fru_a_motif_sum > 0 then flag_fru_a_motif = 1;
            else flag_fru_a_motif  = 0;
            if flag_fru_b_motif_sum > 0 then flag_fru_b_motif = 1;
            else flag_fru_b_motif  = 0;
            if flag_fru_c_motif_sum > 0 then flag_fru_c_motif = 1;
            else flag_fru_c_motif  = 0;
            if flag_fru_a_motif  = 1 or flag_fru_b_motif  = 1 or  flag_fru_c_motif  = 1 then flag_fru_binding = 1;
            else flag_fru_binding = 0;
            if flag_fru_over_sum > 0 then flag_fru_over = 1;
            else flag_fru_over  = 0;
            if flag_dsx_null_sum > 0 then flag_dsx_null = 1;
            else flag_dsx_null  = 0;
            keep primary_fbgn flag_fru_a_motif flag_fru_b_motif
            flag_fru_c_motif flag_fru_binding flag_fru_over flag_dsx_null;
            run;

    proc sort data=ribomerge;
        by primary_fbgn;
        run;

    proc sort data=frufb551;
        by primary_fbgn;
        run;

    data fruMerge;
        merge  ribomerge (in=in1) frufb551 (in=in2);
        by primary_fbgn;
        if in1 and not in2 then do;
            flag_fru_a_motif  = 0;
            flag_fru_b_motif  = 0;
            flag_fru_c_motif  = 0;
            flag_fru_binding  = 0;
            flag_fru_over = 0;
            flag_dsx_null = 0;
        end;
        if in1 then output fruMerge;
        run;

/* Merge on DSX BS */
    libname genelist '!MCLAB/useful_dmel_data/gene_lists/sas_data';

    data dsx;
        set genelist.Luo2011_dsx_binding_site;
        run;

    proc sort data= dsx nodupkey;
        by primary_fbgn;
        run;

    data dsxMerge tocheck;
        merge fruMerge (in=in1) dsx (in=in2);
        by primary_fbgn;
        if in2 then flag_dsx_binding = 1;
        else flag_dsx_binding = 0;
        if in1 and not in2 then flag_damid_sig = 0;
        if in2 and not in1 then output tocheck;
        if in1 then output dsxMerge;
        run;

/* More freqs */
    proc freq data=dsxMerge;
        tables flag_us_fru_dsx*(flag_575 flag_fru_binding flag_dsx_binding) /chisq exact expect;
        run;

/* Export flags */
    proc export data=dsxMerge outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/program_project_genes_us_fru_and_dsx_not_ds_sxl_w_flags.csv' dbms=csv replace;
        putnames=yes;
        run;

/* export Gene list */
    data out;
        set dsxMerge;
        if flag_us_fru_dsx = 1 and flag_fru_over = 0 and flag_fru_binding = 0 and
            flag_dsx_null = 0 and flag_dsx_binding = 0;
        keep primary_fbgn symbol;
        run;


    proc export data=out outfile='!MCLAB/cegs_sem_sd_paper/analysis_output/program_project_gene_list.csv' dbms=csv replace;
        putnames=yes;
        run;



/* Clean up */
    proc datasets ;
        delete ag;
        delete freqs;
        delete fru;
        delete frufb551;
        delete merged;
        delete ribo;
        delete ribomerge;
        delete sums;
        delete tocheck;
        delete us;
        delete us_dsx_fru;
        run; quit;
