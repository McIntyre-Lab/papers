/* 
 * REVISION: 
 *      12/01/2011: Changed FDR level from .05 to .20
 *      12/27/2011: Removed some contrasts
 */

libname fru "!MCLAB/Fru_network/sasdata";
options formdlim='-';

%include "!MCLAB/Fru_network/sas_programs/macro_fdr_mult_plus_indicators_v5.sas";

/* 
 * Because these contrasts labels are so long and I am going to want to
 * transpose these data and use them as column names I decided to simply create
 * a contrast variable that will be used as the column names, then using the
 * idlabel statment I can use the contrast names as labels which can have more
 * then 32 characters
 */

data contrast;
    retain fusion_id contrast source;
    set fru.contrast_by_fusion;
    label probf  = 'p_contrast';
    label fvalue = 'f_contrast';
    rename probf  = p_contrast
           fvalue = f_contrast;
     if source = 'AHCSxAHBerM'                 then contrast =  1 ;
     if source = 'AHCSFemalexAHBerF'           then contrast =  2 ;
     if source = 'AHCSFemalexAHCS'             then contrast =  3 ;
     if source = 'AHBerMxAHBerF'               then contrast =  4 ;
     if source = 'AHCSxAHBerF'                 then contrast =  5 ;
     if source = 'AHCSFemalexAHBerM'           then contrast =  6 ;
     if source = 'AHCSxAHFruP14440'            then contrast =  7 ;
     if source = 'AHBerMxAHFruP14440'          then contrast =  8 ;
     if source = 'AHCSxAHFruW12ChaM5'          then contrast =  9 ;
     if source = 'AHBerMxAHFruW12ChaM5'        then contrast =  10;
     if source = 'AHCSxAHMaleFruM(A)'          then contrast =  11;
     if source = 'AHBerMxAHMaleFruM(A)'        then contrast =  12;
     if source = 'AHCSxAHMaleFruM(B)'          then contrast =  13;
     if source = 'AHBerMxAHMaleFruM(B)'        then contrast =  14;
     if source = 'AHCSxAHMaleFruM(C)'          then contrast =  15;
     if source = 'AHBerMxAHMaleFruM(C)'        then contrast =  16;
     if source = 'AHCSFemalexAHFemaleFruM(A)'  then contrast =  17;
     if source = 'AHBerFxAHFemaleFruM(A)'      then contrast =  18;
     if source = 'AHCSFemalexAHFemaleFruM(B)'  then contrast =  19;
     if source = 'AHBerFxAHFemaleFruM(B)'      then contrast =  20;
     if source = 'AHCSFemalexAHFemaleFruM(C)'  then contrast =  21;
     if source = 'AHBerFxAHFemaleFruM(C)'      then contrast =  22;
     drop dependent;
    run;

proc sort data=contrast;
    by fusion_id contrast;
    run;

/* 
 * Transpose the contrasts each column will have the value p_all_contrast# but
 * will be labled with the contrast label. Then create a flag if the p_value is
 * significant.
 */

proc transpose data=contrast out = p_trans prefix = p_all_ ;
    by fusion_id;
    id contrast;
    var p_contrast;
    run;

data p_trans;
    set p_trans;
    label p_all_1  =  'p_all_AHCSxAHBerM'               ;
    label p_all_2  =  'p_all_AHCSFemalexAHBerF'         ;
    label p_all_3  =  'p_all_AHCSFemalexAHCS'           ;
    label p_all_4  =  'p_all_AHBerMxAHBerF'             ;
    label p_all_5  =  'p_all_AHCSxAHBerF'               ;
    label p_all_6  =  'p_all_AHCSFemalexAHBerM'         ;
    label p_all_7  =  'p_all_AHCSxAHFruP14440'          ;
    label p_all_8  =  'p_all_AHBerMxAHFruP14440'        ;
    label p_all_9  =  'p_all_AHCSxAHFruW12ChaM5'        ;
    label p_all_10 =  'p_all_AHBerMxAHFruW12ChaM5'      ;
    label p_all_11 =  'p_all_AHCSxAHMaleFruM(A)'        ;
    label p_all_12 =  'p_all_AHBerMxAHMaleFruM(A)'      ;
    label p_all_13 =  'p_all_AHCSxAHMaleFruM(B)'        ;
    label p_all_14 =  'p_all_AHBerMxAHMaleFruM(B)'      ;
    label p_all_15 =  'p_all_AHCSxAHMaleFruM(C)'        ;
    label p_all_16 =  'p_all_AHBerMxAHMaleFruM(C)'      ;
    label p_all_17 =  'p_all_AHCSFemalexAHFemaleFruM(A)';
    label p_all_18 =  'p_all_AHBerFxAHFemaleFruM(A)'    ;
    label p_all_19 =  'p_all_AHCSFemalexAHFemaleFruM(B)';
    label p_all_20 =  'p_all_AHBerFxAHFemaleFruM(B)'    ;
    label p_all_21 =  'p_all_AHCSFemalexAHFemaleFruM(C)';
    label p_all_22 =  'p_all_AHBerFxAHFemaleFruM(C)'    ;
    run;

%macro flag_pvals (pval,label);

    data flag_&pval;
        set p_trans;
        if &pval = . then flag_&pval._05 = .;
            else if &pval < 0.05 then flag_&pval._05 = 1;
            else flag_&pval._05 = 0;
        label flag_&pval._05 = "flag_p_all_05_&label.";
        keep fusion_id &pval flag_&pval._05;
        run;
    
    proc sort data = flag_&pval;
        by fusion_id;
        run;

%mend;

proc contents data=p_trans varnum;
run;

%flag_pvals(p_all_1,AHCSxAHBerM);
%flag_pvals(p_all_2,AHCSFemalexAHBerF);
%flag_pvals(p_all_3,AHCSFemalexAHCS);
%flag_pvals(p_all_4,AHBerMxAHBerF);
%flag_pvals(p_all_5,AHCSxAHBerF);
%flag_pvals(p_all_6,AHCSFemalexAHBerM);
%flag_pvals(p_all_7,AHCSxAHFruP14440);
%flag_pvals(p_all_8,AHBerMxAHFruP14440);
%flag_pvals(p_all_9,AHCSxAHFruW12ChaM5);
%flag_pvals(p_all_10,AHBerMxAHFruW12ChaM5);
%flag_pvals(p_all_11,AHCSxAHMaleFruM(A));
%flag_pvals(p_all_12,AHBerMxAHMaleFruM(A));
%flag_pvals(p_all_13,AHCSxAHMaleFruM(B));
%flag_pvals(p_all_14,AHBerMxAHMaleFruM(B));
%flag_pvals(p_all_15,AHCSxAHMaleFruM(C));
%flag_pvals(p_all_16,AHBerMxAHMaleFruM(C));
%flag_pvals(p_all_17,AHCSFemalexAHFemaleFruM(A));
%flag_pvals(p_all_18,AHBerFxAHFemaleFruM(A));
%flag_pvals(p_all_19,AHCSFemalexAHFemaleFruM(B));
%flag_pvals(p_all_20,AHBerFxAHFemaleFruM(B));
%flag_pvals(p_all_21,AHCSFemalexAHFemaleFruM(C));
%flag_pvals(p_all_22,AHBerFxAHFemaleFruM(C));

/*
 * Calculate FDR for the contrasts, then create a flag if the fdr value is
 * significant
 */

proc sort data=contrast ;
    by p_contrast;
    run;

%fdr_mult(contrast,p_contrast, fusion_id, 0.2, 0.1, 0.05);

data fdr_p_contrast;
    set fdr_p_contrast;
    keep fusion_id source contrast fdr_p_contrast;
    label fdr_p_contrast = "fdr_p_contrast";
    run;

proc transpose data=fdr_p_contrast out=fdr_contrast_trans prefix=fdr_p_contrast_;
    by fusion_id;
    id contrast;
    var fdr_p_contrast;
    run;

data fdr_contrast_trans;
    set fdr_contrast_trans;
    label fdr_p_contrast_1  =  'fdr_p_contrast_AHCSxAHBerM'               ;
    label fdr_p_contrast_2  =  'fdr_p_contrast_AHCSFemalexAHBerF'         ;
    label fdr_p_contrast_3  =  'fdr_p_contrast_AHCSFemalexAHCS'           ;
    label fdr_p_contrast_4  =  'fdr_p_contrast_AHBerMxAHBerF'             ;
    label fdr_p_contrast_5  =  'fdr_p_contrast_AHCSxAHBerF'               ;
    label fdr_p_contrast_6  =  'fdr_p_contrast_AHCSFemalexAHBerM'         ;
    label fdr_p_contrast_7  =  'fdr_p_contrast_AHCSxAHFruP14440'          ;
    label fdr_p_contrast_8  =  'fdr_p_contrast_AHBerMxAHFruP14440'        ;
    label fdr_p_contrast_9  =  'fdr_p_contrast_AHCSxAHFruW12ChaM5'        ;
    label fdr_p_contrast_10 =  'fdr_p_contrast_AHBerMxAHFruW12ChaM5'      ;
    label fdr_p_contrast_11 =  'fdr_p_contrast_AHCSxAHMaleFruM(A)'        ;
    label fdr_p_contrast_12 =  'fdr_p_contrast_AHBerMxAHMaleFruM(A)'      ;
    label fdr_p_contrast_13 =  'fdr_p_contrast_AHCSxAHMaleFruM(B)'        ;
    label fdr_p_contrast_14 =  'fdr_p_contrast_AHBerMxAHMaleFruM(B)'      ;
    label fdr_p_contrast_15 =  'fdr_p_contrast_AHCSxAHMaleFruM(C)'        ;
    label fdr_p_contrast_16 =  'fdr_p_contrast_AHBerMxAHMaleFruM(C)'      ;
    label fdr_p_contrast_17 =  'fdr_p_contrast_AHCSFemalexAHFemaleFruM(A)';
    label fdr_p_contrast_18 =  'fdr_p_contrast_AHBerFxAHFemaleFruM(A)'    ;
    label fdr_p_contrast_19 =  'fdr_p_contrast_AHCSFemalexAHFemaleFruM(B)';
    label fdr_p_contrast_20 =  'fdr_p_contrast_AHBerFxAHFemaleFruM(B)'    ;
    label fdr_p_contrast_21 =  'fdr_p_contrast_AHCSFemalexAHFemaleFruM(C)';
    label fdr_p_contrast_22 =  'fdr_p_contrast_AHBerFxAHFemaleFruM(C)'    ;
    run;


%macro flag_fdrs (fdr,label);

    data flag_&fdr;
        set fdr_contrast_trans;
        if &fdr = . then flag_&fdr._20= .;
            else if &fdr < 0.20 then flag_&fdr._20 = 1;
            else flag_&fdr._20 = 0;
        label flag_&fdr._20 = "flag_fdr_20_&label.";
        keep fusion_id &fdr flag_&fdr._20;
        run;

    proc sort data = flag_&fdr;
        by fusion_id;
        run;

%mend ;

proc contents data=fdr_contrast_trans varnum;
run;

%flag_fdrs(fdr_p_contrast_1,AHCSxAHFruP14440);
%flag_fdrs(fdr_p_contrast_1,AHCSxAHBerM);
%flag_fdrs(fdr_p_contrast_2,AHCSFemalexAHBerF);
%flag_fdrs(fdr_p_contrast_3,AHCSFemalexAHCS);
%flag_fdrs(fdr_p_contrast_4,AHBerMxAHBerF);
%flag_fdrs(fdr_p_contrast_5,AHCSxAHBerF);
%flag_fdrs(fdr_p_contrast_6,AHCSFemalexAHBerM);
%flag_fdrs(fdr_p_contrast_7,AHCSxAHFruP14440);
%flag_fdrs(fdr_p_contrast_8,AHBerMxAHFruP14440);
%flag_fdrs(fdr_p_contrast_9,AHCSxAHFruW12ChaM5);
%flag_fdrs(fdr_p_contrast_10,AHBerMxAHFruW12ChaM5);
%flag_fdrs(fdr_p_contrast_11,AHCSxAHMaleFruM(A));
%flag_fdrs(fdr_p_contrast_12,AHBerMxAHMaleFruM(A));
%flag_fdrs(fdr_p_contrast_13,AHCSxAHMaleFruM(B));
%flag_fdrs(fdr_p_contrast_14,AHBerMxAHMaleFruM(B));
%flag_fdrs(fdr_p_contrast_15,AHCSxAHMaleFruM(C));
%flag_fdrs(fdr_p_contrast_16,AHBerMxAHMaleFruM(C));
%flag_fdrs(fdr_p_contrast_17,AHCSFemalexAHFemaleFruM(A));
%flag_fdrs(fdr_p_contrast_18,AHBerFxAHFemaleFruM(A));
%flag_fdrs(fdr_p_contrast_19,AHCSFemalexAHFemaleFruM(B));
%flag_fdrs(fdr_p_contrast_20,AHBerFxAHFemaleFruM(B));
%flag_fdrs(fdr_p_contrast_21,AHCSFemalexAHFemaleFruM(C));
%flag_fdrs(fdr_p_contrast_22,AHBerFxAHFemaleFruM(C));

/*
 * Merge everything together 
 */

data fru.flag_fdr_contrast_by_fusion;
    merge flag_p_all_1 flag_p_all_2 flag_p_all_3 flag_p_all_4 
          flag_p_all_5 flag_p_all_6 flag_p_all_7 flag_p_all_8 
          flag_p_all_9 flag_p_all_10 flag_p_all_11 flag_p_all_12 
          flag_p_all_13 flag_p_all_14 flag_p_all_15 flag_p_all_16 
          flag_p_all_17 flag_p_all_18 flag_p_all_19 flag_p_all_20 
          flag_p_all_21 flag_p_all_22

          flag_fdr_p_contrast_1 flag_fdr_p_contrast_2 flag_fdr_p_contrast_3 
          flag_fdr_p_contrast_4 flag_fdr_p_contrast_5 flag_fdr_p_contrast_6 
          flag_fdr_p_contrast_7 flag_fdr_p_contrast_8 flag_fdr_p_contrast_9 
          flag_fdr_p_contrast_10 flag_fdr_p_contrast_11 flag_fdr_p_contrast_12 
          flag_fdr_p_contrast_13 flag_fdr_p_contrast_14 flag_fdr_p_contrast_15 
          flag_fdr_p_contrast_16 flag_fdr_p_contrast_17 flag_fdr_p_contrast_18 
          flag_fdr_p_contrast_19 flag_fdr_p_contrast_20 flag_fdr_p_contrast_21 
          flag_fdr_p_contrast_22;
    by fusion_id;
    run;
