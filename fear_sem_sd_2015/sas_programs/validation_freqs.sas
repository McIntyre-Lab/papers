/*******************************************************************************
* Filename: validation_freqs.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: 
*
*******************************************************************************/

/* Libraries
libname sem '!MCLAB/cegs_sem_sd_paper/sas_data';
libname genelist '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_sem_sd_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

data valid;
    set SEM.validation_set;
    run;

data valid2;
    set valid;
    drop model modelnum bic path;
    run;

proc sort data=valid2 nodupkey;
    by primary_fbgn;
    run;

proc sort data=valid2 nodupkey;
    by primary_fbgn;
    run;

proc freq data=valid2;
    tables flag_grn_expansion;
    tables flag_grn_expansion*(flag_dsxNullf_repressed
    flag_chang_ds_tra
    flag_chang_ups_tra
    flag_chang_tra_bs
    flag_goldman_ds_tra);
    run;

data check;
    set valid2;
    if flag_chang_tra_bs=1 and flag_grn_expansion=1;
    run;

proc freq data=check;
    tables flag_grn_expansion*(flag_dsxNullf_repressed
    flag_chang_ds_tra
    flag_chang_ups_tra
    flag_chang_tra_bs
    flag_goldman_ds_tra);
    run;

/* Pull All models and figure out what top one are */
%macro pull_models(name, gene);
    data &name;
        set valid;
        where symbol ? "&gene";
        run;

    proc sort data=&name;
        by bic;
        run;
%mend;
%pull_models(Psa ,Psa ) ;
%pull_models(msl2 ,msl-2 ) ;
%pull_models(Mob ,Mob ) ; *yes;
%pull_models(B52 ,B52 ) ; *yes;
%pull_models(sqd ,sqd ) ; *yes;
%pull_models(Psi ,Psi ) ;*no;
%pull_models(mub ,mub ) ; *not in merged?;
%pull_models(Rm62 ,Rm62 ) ; *winner;
%pull_models(ps ,ps ) ; *no evidence;
%pull_models(hts ,hts ) ; *downstream sxl;
%pull_models(ltd ,ltd ) ; *wahoo;
%pull_models(garz ,garz ) ; *no;
%pull_models(jdp ,jdp ) ; *no;
%pull_models(bbc ,bbc ) ; *no;
%pull_models(CG4662, CG4662 ) ; *no;

/* Look at Freqs */
    proc freq data=valid2;
        tables flag_chang_ds_tra*flag_chang_ups_tra;
        run; * 635 genes upstream, 350 genes downstream of tra;

    proc freq data=valid2;
        where flag_chang_ds_tra=1 or flag_chang_ups_tra=1;
        tables flag_chang_tra_bs;
        run; * 2 genes contained tra bs (of 985);

    proc freq data=valid2;
        where flag_chang_ds_tra=1 or flag_chang_ups_tra=1;
        tables flag_goldman_ds_tra;
        run; * 39 downstream in goldman, 3 not downstream in goldman, 943 not in goldman;

    data check2;
        set valid2;
        if flag_goldman_ds_tra=1 and (flag_chang_ds_tra=1 or flag_chang_ups_tra=1);
        run; *39 obs, ok;




data splice;
    set valid;
    if psex_by_probe ne .;
    run;

proc sort data=splice;
    by symbol bic;
    run;

proc freq data=splice;
    tables symbol/out=count_gene;
    run;

data check_best;
    set splice;
    by symbol;
    if first.symbol;
    run;

proc freq data=check_best;
    tables model;
    run;

/*
Model_1                  7       15.22             7
Model_10                 1        2.17             8
Model_22                 1        2.17             9
Model_24                 1        2.17            10
Model_29                 2        4.35            12
Model_31                 1        2.17            13
Model_32                 2        4.35            15
Model_35                 3        6.52            18
Model_36                 2        4.35            20
Model_Baseline          26       56.52            46

20 with model / 26 select baseline 

*/




/*
Fear 2014:
flag_sem_added_gene = 1 if the SEM model added the gene to the SD pathway
flag_dsxNull_induced = 1 if a gene had increased expression when dsx was knocked out
flag_dsxNull_repressed= 1 if a gene had decreased expression when dsx was knocked out (since dsx is a TF I could guess these are more likely direct targets)
 
Chang 2011:
flag_chang_female_bias = 1 if there was a diff between MvsF and it was biased toward females
flag_chang_male_bias= 1 if there was a diff between MvsF and it was biased toward males
flag_chang_tra_bias = 1 if there was a diff between FvsTra pseudo males and it was biased toward tra pseudomales
flag_chang_traf_bias = 1 if there was a diff between FvsTra pseudo males and it was biased toward females
flag_chang_ds_tra = 1 if looking at the results the gene appears to be downstream of tra
flag_chang_ups_tra = 1 if looking at the results the gene appears to be upstream of tra
flag_tra_bs = 1 if contains a tra dna binding site
 
Goldman 2007:
flag_goldman_female_bias = 1 if there was a diff between MvsF and it was biased toward females
flag_goldman_male_bias bias = 1 if there was a diff between MvsF and it was biased toward males
flag_goldman_ds_tra =  1 if looking at the results the gene appears to be downstream of tra
*/
