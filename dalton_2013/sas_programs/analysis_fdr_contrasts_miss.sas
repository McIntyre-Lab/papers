/* 
 * REVISION: 
 *      12/01/2011: Changed FDR level from .05 to .20
 *      12/27/2011: Removed some contrasts
 */

options formdlim='-';

%include "!MCLAB/arbeitman/arbeitman_fru_network/sas_programs/macro_fdr_mult_plus_indicators_v5.sas";

/* 
 * Because these contrasts labels are so long and I am going to want to
 * transpose these data and use them as column names I decided to simply create
 * a contrast variable that will be used as the column names, then using the
 * idlabel statment I can use the contrast names as labels which can have more
 * then 32 characters
 */

data contrast;
    retain fusion_id contrast source;
    set fru.contrast_by_fusion_miss;
    label probf  = 'p_miss';
    label fvalue = 'f_miss';
    rename probf  = p_miss
           fvalue = f_miss;
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

proc transpose data=contrast out = p_trans prefix = p_miss_ ;
    by fusion_id;
    id contrast;
    var p_miss;
    run;

data p_trans;
    set p_trans;
    label p_miss_1  =  'p_miss_AHCSxAHBerM'               ;
    label p_miss_2  =  'p_miss_AHCSFemalexAHBerF'         ;
    label p_miss_3  =  'p_miss_AHCSFemalexAHCS'           ;
    label p_miss_4  =  'p_miss_AHBerMxAHBerF'             ;
    label p_miss_5  =  'p_miss_AHCSxAHBerF'               ;
    label p_miss_6  =  'p_miss_AHCSFemalexAHBerM'         ;
    label p_miss_7  =  'p_miss_AHCSxAHFruP14440'          ;
    label p_miss_8  =  'p_miss_AHBerMxAHFruP14440'        ;
    label p_miss_9  =  'p_miss_AHCSxAHFruW12ChaM5'        ;
    label p_miss_10 =  'p_miss_AHBerMxAHFruW12ChaM5'      ;
    label p_miss_11 =  'p_miss_AHCSxAHMaleFruM(A)'        ;
    label p_miss_12 =  'p_miss_AHBerMxAHMaleFruM(A)'      ;
    label p_miss_13 =  'p_miss_AHCSxAHMaleFruM(B)'        ;
    label p_miss_14 =  'p_miss_AHBerMxAHMaleFruM(B)'      ;
    label p_miss_15 =  'p_miss_AHCSxAHMaleFruM(C)'        ;
    label p_miss_16 =  'p_miss_AHBerMxAHMaleFruM(C)'      ;
    label p_miss_17 =  'p_miss_AHCSFemalexAHFemaleFruM(A)';
    label p_miss_18 =  'p_miss_AHBerFxAHFemaleFruM(A)'    ;
    label p_miss_19 =  'p_miss_AHCSFemalexAHFemaleFruM(B)';
    label p_miss_20 =  'p_miss_AHBerFxAHFemaleFruM(B)'    ;
    label p_miss_21 =  'p_miss_AHCSFemalexAHFemaleFruM(C)';
    label p_miss_22 =  'p_miss_AHBerFxAHFemaleFruM(C)'    ;
    run;

%macro flag_pvals (pval,label);

    data flag_&pval;
        set p_trans;
        if &pval = . then flag_&pval._05 = .;
            else if &pval < 0.05 then flag_&pval._05 = 1;
            else flag_&pval._05 = 0;
        label flag_&pval._05 = "flag_p_miss_05_&label.";
        keep fusion_id &pval flag_&pval._05;
        run;
    
    proc sort data = flag_&pval;
        by fusion_id;
        run;

%mend;

proc contents data=p_trans varnum;
run;

%flag_pvals(p_miss_1,AHCSxAHBerM);
%flag_pvals(p_miss_2,AHCSFemalexAHBerF);
%flag_pvals(p_miss_3,AHCSFemalexAHCS);
%flag_pvals(p_miss_4,AHBerMxAHBerF);
%flag_pvals(p_miss_5,AHCSxAHBerF);
%flag_pvals(p_miss_6,AHCSFemalexAHBerM);
%flag_pvals(p_miss_7,AHCSxAHFruP14440);
%flag_pvals(p_miss_8,AHBerMxAHFruP14440);
%flag_pvals(p_miss_9,AHCSxAHFruW12ChaM5);
%flag_pvals(p_miss_10,AHBerMxAHFruW12ChaM5);
%flag_pvals(p_miss_11,AHCSxAHMaleFruM(A));
%flag_pvals(p_miss_12,AHBerMxAHMaleFruM(A));
%flag_pvals(p_miss_13,AHCSxAHMaleFruM(B));
%flag_pvals(p_miss_14,AHBerMxAHMaleFruM(B));
%flag_pvals(p_miss_15,AHCSxAHMaleFruM(C));
%flag_pvals(p_miss_16,AHBerMxAHMaleFruM(C));
%flag_pvals(p_miss_17,AHCSFemalexAHFemaleFruM(A));
%flag_pvals(p_miss_18,AHBerFxAHFemaleFruM(A));
%flag_pvals(p_miss_19,AHCSFemalexAHFemaleFruM(B));
%flag_pvals(p_miss_20,AHBerFxAHFemaleFruM(B));
%flag_pvals(p_miss_21,AHCSFemalexAHFemaleFruM(C));
%flag_pvals(p_miss_22,AHBerFxAHFemaleFruM(C));

/*
 * Calculate FDR for the contrasts, then create a flag if the fdr value is
 * significant
 */

proc sort data=contrast ;
    by p_contrast;
    run;

%fdr_mult(contrast,p_miss, fusion_id, 0.2, 0.1, 0.05);

data fdr_p_miss;
    set fdr_p_miss;
    keep fusion_id source contrast fdr_p_miss;
    label fdr_p_miss = "fdr_p_miss";
    run;

proc transpose data=fdr_p_miss out=fdr_miss_trans prefix=fdr_p_miss_;
    by fusion_id;
    id contrast;
    var fdr_p_miss;
    run;

data fdr_miss_trans;
    set fdr_miss_trans;
    label fdr_p_miss_1  =  'fdr_p_miss_AHCSxAHBerM'               ;
    label fdr_p_miss_2  =  'fdr_p_miss_AHCSFemalexAHBerF'         ;
    label fdr_p_miss_3  =  'fdr_p_miss_AHCSFemalexAHCS'           ;
    label fdr_p_miss_4  =  'fdr_p_miss_AHBerMxAHBerF'             ;
    label fdr_p_miss_5  =  'fdr_p_miss_AHCSxAHBerF'               ;
    label fdr_p_miss_6  =  'fdr_p_miss_AHCSFemalexAHBerM'         ;
    label fdr_p_miss_7  =  'fdr_p_miss_AHCSxAHFruP14440'          ;
    label fdr_p_miss_8  =  'fdr_p_miss_AHBerMxAHFruP14440'        ;
    label fdr_p_miss_9  =  'fdr_p_miss_AHCSxAHFruW12ChaM5'        ;
    label fdr_p_miss_10 =  'fdr_p_miss_AHBerMxAHFruW12ChaM5'      ;
    label fdr_p_miss_11 =  'fdr_p_miss_AHCSxAHMaleFruM(A)'        ;
    label fdr_p_miss_12 =  'fdr_p_miss_AHBerMxAHMaleFruM(A)'      ;
    label fdr_p_miss_13 =  'fdr_p_miss_AHCSxAHMaleFruM(B)'        ;
    label fdr_p_miss_14 =  'fdr_p_miss_AHBerMxAHMaleFruM(B)'      ;
    label fdr_p_miss_15 =  'fdr_p_miss_AHCSxAHMaleFruM(C)'        ;
    label fdr_p_miss_16 =  'fdr_p_miss_AHBerMxAHMaleFruM(C)'      ;
    label fdr_p_miss_17 =  'fdr_p_miss_AHCSFemalexAHFemaleFruM(A)';
    label fdr_p_miss_18 =  'fdr_p_miss_AHBerFxAHFemaleFruM(A)'    ;
    label fdr_p_miss_19 =  'fdr_p_miss_AHCSFemalexAHFemaleFruM(B)';
    label fdr_p_miss_20 =  'fdr_p_miss_AHBerFxAHFemaleFruM(B)'    ;
    label fdr_p_miss_21 =  'fdr_p_miss_AHCSFemalexAHFemaleFruM(C)';
    label fdr_p_miss_22 =  'fdr_p_miss_AHBerFxAHFemaleFruM(C)'    ;
    run;


%macro flag_fdrs (fdr,label);

    data flag_&fdr;
        set fdr_miss_trans;
        if &fdr = . then flag_&fdr._20= .;
            else if &fdr < 0.20 then flag_&fdr._20 = 1;
            else flag_&fdr._20 = 0;
        label flag_&fdr._tmp = "flag_fdr_miss_&label.";
        keep fusion_id &fdr flag_&fdr._20;
        run;

    proc sort data = flag_&fdr;
        by fusion_id;
        run;

%mend ;

proc contents data=fdr_contrast_trans varnum;
run;

%flag_fdrs(fdr_p_miss_1,AHCSxAHFruP14440);
%flag_fdrs(fdr_p_miss_1,AHCSxAHBerM);
%flag_fdrs(fdr_p_miss_2,AHCSFemalexAHBerF);
%flag_fdrs(fdr_p_miss_3,AHCSFemalexAHCS);
%flag_fdrs(fdr_p_miss_4,AHBerMxAHBerF);
%flag_fdrs(fdr_p_miss_5,AHCSxAHBerF);
%flag_fdrs(fdr_p_miss_6,AHCSFemalexAHBerM);
%flag_fdrs(fdr_p_miss_7,AHCSxAHFruP14440);
%flag_fdrs(fdr_p_miss_8,AHBerMxAHFruP14440);
%flag_fdrs(fdr_p_miss_9,AHCSxAHFruW12ChaM5);
%flag_fdrs(fdr_p_miss_10,AHBerMxAHFruW12ChaM5);
%flag_fdrs(fdr_p_miss_11,AHCSxAHMaleFruM(A));
%flag_fdrs(fdr_p_miss_12,AHBerMxAHMaleFruM(A));
%flag_fdrs(fdr_p_miss_13,AHCSxAHMaleFruM(B));
%flag_fdrs(fdr_p_miss_14,AHBerMxAHMaleFruM(B));
%flag_fdrs(fdr_p_miss_15,AHCSxAHMaleFruM(C));
%flag_fdrs(fdr_p_miss_16,AHBerMxAHMaleFruM(C));
%flag_fdrs(fdr_p_miss_17,AHCSFemalexAHFemaleFruM(A));
%flag_fdrs(fdr_p_miss_18,AHBerFxAHFemaleFruM(A));
%flag_fdrs(fdr_p_miss_19,AHCSFemalexAHFemaleFruM(B));
%flag_fdrs(fdr_p_miss_20,AHBerFxAHFemaleFruM(B));
%flag_fdrs(fdr_p_miss_21,AHCSFemalexAHFemaleFruM(C));
%flag_fdrs(fdr_p_miss_22,AHBerFxAHFemaleFruM(C));

/*
 * Merge everything together 
 */

data fru.flag_fdr_contrast_by_fusion_miss;
    merge flag_p_miss_1 flag_p_miss_2 flag_p_miss_3 flag_p_miss_4 
          flag_p_miss_5 flag_p_miss_6 flag_p_miss_7 flag_p_miss_8 
          flag_p_miss_9 flag_p_miss_10 flag_p_miss_11 flag_p_miss_12 
          flag_p_miss_13 flag_p_miss_14 flag_p_miss_15 flag_p_miss_16 
          flag_p_miss_17 flag_p_miss_18 flag_p_miss_19 flag_p_miss_20 
          flag_p_miss_21 flag_p_miss_22

          flag_fdr_p_miss_1 flag_fdr_p_miss_2 flag_fdr_p_miss_3 
          flag_fdr_p_miss_4 flag_fdr_p_miss_5 flag_fdr_p_miss_6 
          flag_fdr_p_miss_7 flag_fdr_p_miss_8 flag_fdr_p_miss_9 
          flag_fdr_p_miss_10 flag_fdr_p_miss_11 flag_fdr_p_miss_12 
          flag_fdr_p_miss_13 flag_fdr_p_miss_14 flag_fdr_p_miss_15 
          flag_fdr_p_miss_16 flag_fdr_p_miss_17 flag_fdr_p_miss_18 
          flag_fdr_p_miss_19 flag_fdr_p_miss_20 flag_fdr_p_miss_21 
          flag_fdr_p_miss_22;
    by fusion_id;
    run;
