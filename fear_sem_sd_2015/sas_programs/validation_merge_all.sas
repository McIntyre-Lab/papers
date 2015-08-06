/*******************************************************************************
* Filename: validation_merge_all.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: There are several lines of external data that we are using for
* validation. Merge all of these datasets onto the 1390 genes added by the SEM.
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname dsx '!MCLAB/arbeitman/arbeitman_dsx/sas_data';
libname fru '!MCLAB/arbeitman/arbeitman_fru_network/sasdata';
libname dmel551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';
libname genelist '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge BICs to FBgn */
    proc sort data=SEM.cegsV_gene2fbgn;
        by gene;
        run; * 8822 obs;

    proc sort data=SEM.cegsV_ag_yp2_stack_bic;
        by gene;
        run; * 334780 obs;

    proc sort data=SEM.cegsV_ag_yp2_added_genes;
        by gene;
        run; * 1390 obs;

    data agmerge;
        merge SEM.cegsV_gene2fbgn (in=in1) SEM.cegsV_ag_yp2_stack_bic (in=in2) 
              SEM.cegsV_ag_yp2_added_genes (in=in3);
        by gene;
        if in3 then flag_grn_expansion = 1;
        else flag_grn_expansion = 0;
        if in2;
        drop gene;
        run; * 334780 obs;

    /* Check and make sure that I have 1390 genes with added path

    data mycheck;
        set agmerge;
        where flag_grn_expansion = 1;
        run;

    proc sort data=mycheck nodupkey;
        by primary_fbgn;
        run; 1390 obs;

    */

/* Merge to Model Design */
    proc sort data=agmerge;
        by model;
        run;

    proc sort data=SEM.cegsV_ag_model_design_file;
        by model;
        run;

    data design;
        merge agmerge (in=in1) SEM.cegsV_ag_model_design_file (in=in2);
        by model;
        if in1;
        run;

/* dsxNUll gene lists to fb5.51 to check fbgns */
    * Merge to fb530;
    data fbgn;
        set DMEL530.symbol2fbgn;
        run; * 228781 obs;

    data dsx_ctrl_female_induced;
        length fbgn_cat $11.;
        set SEM.dsx_ctrl_female_induced;
        rename fbgn_cat = primary_fbgn;
        run; * 246 obs;

    data dsx_ctrl_female_repressed;
        length fbgn_cat $11.;
        set SEM.dsx_ctrl_female_repressed;
        rename fbgn_cat = primary_fbgn;
        run; * 205 obs;

    data dsxNullf_induced;
        length fbgn_cat $11.;
        set SEM.dsxNullf_induced;
        rename fbgn_cat = primary_fbgn;
        run; * 340 obs;

    data dsxNullf_repressed;
        length fbgn_cat $11.;
        set SEM.dsxNullf_repressed;
        rename fbgn_cat = primary_fbgn;
        run; * 208 obs;

    * Also add the fru binding site flags;
    data flag_dalton_fru_bs;
        set FRU.Motif_flags_and_cnts;
        sums =  flag_fru_a_motif + flag_fru_b_motif + flag_fru_c_motif;
        if sums > 0 then flag_dalton_fru_bs = 1;
        else flag_dalton_fru_bs = 0;
        keep primary_fbgn flag_dalton_fru_bs;
        run; * 14903 obs;
    
    proc sort data=fbgn;
        by primary_fbgn;
        run;

    proc sort data=dsx_ctrl_female_induced;
        by primary_fbgn;
        run;

    proc sort data=dsx_ctrl_female_repressed;
        by primary_fbgn;
        run;

    proc sort data=dsxNullf_induced;
        by primary_fbgn;
        run;

    proc sort data=dsxNullf_repressed;
        by primary_fbgn;
        run;

     proc sort data=flag_dalton_fru_bs;
         by primary_fbgn;
         run;

    data anno oops;
        merge fbgn (in=in1) dsx_ctrl_female_induced (in=in2) dsx_ctrl_female_repressed (in=in3)
              dsxNullf_induced (in=in4) dsxNullf_repressed (in=in5) flag_dalton_fru_bs (in=in6);
        by primary_fbgn;
        if in2 then flag_dsxCtrl_female_induced = 1; else flag_dsxCtrl_female_induced = 0;
        if in3 then flag_dsxCtrl_female_repressed = 1; else flag_dsxCtrl_female_repressed = 0;
        if in4 then flag_dsxNullf_induced = 1; else flag_dsxNullf_induced = 0;
        if in5 then flag_dsxNullf_repressed = 1; else flag_dsxNullf_repressed = 0;
        * Had some merge issues below, some genes changed symbols between fb5.30
        * and fb5.51, used flymine to figure this out;
        if symbol eq 'CG42256' then symbol = 'Dscam2';
        if symbol eq 'Spn1' then symbol = 'SPN1';
        if symbol eq 'Gap1' then symbol = 'GAP1';
        if symbol eq 'CG6199' then symbol = 'Plod '; 
        if symbol eq 'Anxb11' then symbol = 'AnxB11';
        if symbol eq 'CG11299' then symbol = 'Sesn '; 
        if symbol eq 'CG14375' then symbol = 'CCHa2 ';
        if symbol eq 'CG7077' then symbol = 'rdhB ';
        if symbol eq 'Tbp-1' then symbol = 'tbp-1 ';
        if symbol eq 'CG9542' then symbol = 'KFase ';
        if symbol eq 'CG6783' then symbol = 'fabp ';
        if symbol eq 'CG15828' then symbol = 'Apoltp ';
        if symbol eq 'flower' then symbol = 'fwe ';
        if symbol eq 'CG4688' then symbol = 'GstE14 ';
        if symbol eq 'CG6859' then symbol = 'Pex3 ';
        if symbol eq 'CG3244' then symbol = 'Clect27 ';
        if symbol eq 'CG31665' then symbol = 'wry ';
        if symbol eq 'elk' then symbol = 'Elk ';
        if symbol eq 'CG10267' then symbol = 'Zif ';
        if symbol eq 'CG30115' then symbol = 'GEFmeso ';
        if symbol eq 'CG9027' then symbol = 'Sod3 ';
        if symbol eq 'CG31150' then symbol = 'cv-d ';
        if symbol eq 'CG6461' then symbol = 'Ggt-1 ';
        if symbol eq 'CG1908' then symbol = 'Lint-1 ';
        if symbol eq 'CG10990' then symbol = 'Pdcd4 ';
        if symbol eq 'CG6361' then symbol = 'Hayan ';
        if symbol eq 'CG3759' then symbol = 'Mco1 ';
        if symbol eq 'pdfr' then symbol = 'Pdfr ';
        if symbol eq 'CG3056' then symbol = 'ssx ';
        if symbol eq 'CG9619' then symbol = 'Gbs-76A ';
        if symbol eq 'yu' then symbol = 'Yu ';
        if symbol eq 'cals' then symbol = 'Cals ';
        if symbol eq 'CG6719' then symbol = 'mgr ';
        if symbol eq 'Dh' then symbol = 'Dh ';
        if symbol eq 'nimC1' then symbol = 'NimC1 ';
        if symbol eq 'CG6776' then symbol = 'GstO3 ';
        if symbol eq 'Pros54' then symbol = 'Rpn10 ';
        if symbol eq 'l(2)05510' then symbol = 'qsm ';
        if symbol eq 'CG6043' then symbol = 'CG44085 ';
        if symbol eq 'Pros29' then symbol = 'Prosalpha3 ';
        if symbol eq 'CG34400' then symbol = 'dysc ';
        if symbol eq 'G-salpha60A' then symbol = 'Galphas ';
        if symbol eq 'AnnIX' then symbol = 'AnxB9 ';
        if symbol eq 'abba' then symbol = 'tn ';
        if symbol eq 'CG3563' then symbol = 'jvl ';
        if symbol eq 'CG31714' then symbol = 'CG44153 ';
        if symbol eq 'CR11700' then symbol = 'CG11700 ';
        if symbol eq 'CG15928' then symbol = 'cac ';
        if symbol eq 'Pros35' then symbol = 'Prosalpha6 ';
        if symbol eq 'SPN1' then symbol = 'Spn42Dd ';
        if symbol eq 'Mov34' then symbol = 'Rpn8 ';
        if symbol eq 'dro4' then symbol = 'Drsl4 ';
        if symbol eq 'tbp-1' then symbol = 'Rpt5 ';
        if symbol eq 'GAP1' then symbol = 'RasGAP1 ';
        if symbol eq 'dro5' then symbol = 'Drsl5 ';
        if symbol eq 'TepIV' then symbol = 'Tep4 ';
        if symbol eq 'Smg1' then symbol = 'nonC ';
        if symbol eq 'Pros28.1' then symbol = 'Prosalpha4 ';
        if symbol eq 'AnnX' then symbol = 'AnxB10 ';
        if symbol eq 'TepII' then symbol = 'Tep2 ';
        if symbol eq 'nAcRbeta-64B' then symbol = 'nAcRbeta-64B';
        if symbol eq 'Dh' then symbol = 'Dh44 ';
        if in1 then output anno;  * 228781 obs;
        else output oops; * 0 obs;
        drop primary_fbgn;
        run;

    * Merge to fb551;
    data fb551;
        length symbol $56.;
        set DMEL551.symbol2fbgn;
        run; * 229317 obs;

    proc sort data=anno;
        by symbol;
        run; * 228781 obs;

    proc sort data=fb551;
        by symbol;
        run; * 229317 obs;

    data anno2;
        merge fb551 (in=in1) anno (in=in2);
        by symbol;
        if in2 then output anno2;
        drop symbol;
        run; * 228781 obs;

/* Merge Dsx to model */
    proc sort data=anno2;
        by primary_fbgn;
        run;

    proc sort data=design;
        by primary_fbgn;
        run;

    data dsxNull;
        merge design (in=in1) anno2 (in=in2);
        by primary_fbgn;
        if in1 and not in2 then do;
            flag_dsxCtrl_female_induced = 0;
            flag_dsxCtrl_female_repressed = 0;
            flag_dsxNullf_induced = 0;
            flag_dsxNullf_repressed = 0;
        end;
        if in1 then output;
        run; * 334780 obs;

/* Merge on Luo dsx binding sites */
    proc sort data=GENELIST.Luo2011_dsx_binding_site;
        by primary_fbgn;
        run; * 56 obs;

    data luo;
        merge dsxNull (in=in1) GENELIST.Luo2011_dsx_binding_site (in=in2);
        by primary_fbgn;
        if in2 then flag_luo_binding_site = 1;
        else do;
            flag_luo_binding_site = 0;
            flag_damid_sig = 0;
        end;
        drop symbol;
        run; * 334795 obs;

/* Merge on Chang */
    proc sort data=GENELIST.Chang2011_tra;
        by primary_fbgn;
        run; * 8896 obs;

    data chang oops;
        merge luo (in=in1) GENELIST.Chang2011_tra (in=in2);
        by primary_fbgn;
        if in1 then output chang; * 334795 obs;
        else output oops; * 1137 obs;
        run;

/* Merge on McIntyre */
    proc sort data=GENELIST.McIntyre2006;
        by primary_fbgn;
        run; * 305 obs;

    data mcintyre oops;
        merge chang (in=in1) GENELIST.McIntyre2006 (in=in2);
        by primary_fbgn;
        if in1 then output mcintyre; * 334795 obs;
        else output oops; * 73 obs;
        run;

    data mcintyre;
        set mcintyre;
        if psex_by_probe eq . then flag_mcintyre_sex_bias_splice = .;
        else if psex_by_probe le .05 then flag_mcintyre_sex_bias_splice = 1;
        else flag_mcintyre_sex_bias_splice = 0;
        run;

/* Merge on Goldman */
    proc sort data=GENELIST.Goldman2007_tra;
        by primary_fbgn;
        run; * 135 obs;

    data goldman oops;
        merge mcintyre (in=in1) GENELIST.Goldman2007_tra (in=in2);
        by primary_fbgn;
        if in1 then output goldman; * 334795 obs;
        else output oops; * 36 obs;
        run;

/* Merge on gene symbol */
    proc sort data = DMEL551.symbol2fbgn;
        by primary_fbgn;
        run; * 229317 obs;

    proc sort data=goldman;
        by primary_fbgn;
        run;

    data merged oops;
        retain primary_fbgn symbol;
        merge goldman (in=in1) DMEL551.symbol2fbgn (in=in2);
        by primary_fbgn;
        if in1 then output merged; * 334795 obs;
        else output oops; * 220492 obs;
        run;
        
/* Make Perm Dataset */
    data SEM.validation_set;
        set merged;
        run; * 334795 obs;

/* Clean Up */
    proc datasets nolist;
        delete agmerge;
        delete anno;
        delete anno2;
        delete chang;
        delete design;
        delete dsxnull;
        delete dsxnullf_induced;
        delete dsxnullf_repressed;
        delete dsx_ctrl_female_induced;
        delete dsx_ctrl_female_repressed;
        delete fb551;
        delete fbgn;
        delete goldman;
        delete luo;
        delete mcintyre;
        delete oops;
        delete merged;
        run; quit;
