/*
 * REVISIONS: 12/23/2011
 *             - Removed all Results information and moved to NOTEBOOK
 *             - Changed Contrasts to only include contrasts for Adult Head
 */
libname fru "!MCLAB/Fru_network/sasdata";

proc sort data=fru.foranalysis_no_multi;
    by symbol_cat;
    run;

proc freq data=fru.foranalysis_no_multi;
    tables trt;
    run;

/*
  1   AH_BerF
  2   AH_BerM
  3   AH_CS
  4   AH_CSFemale
  5   AH_Female_FruM(A)
  6   AH_Female_FruM(B)
  7   AH_Female_FruM(C)
  8   AH_FruP14_440
  9   AH_FruW12_ChaM5
  10  AH_Male_FruM(A)
  11  AH_Male_FruM(B)
  12  AH_Male_FruM(C)
*/

ods listing close;

proc glm data=fru.foranalysis_no_multi plots=none;
    by symbol_cat;
    class trt fusion_id;
    model logrpkm = fusion_id trt fusion_id*trt ;
    lsmeans fusion_id trt fusion_id fusion_id*trt/pdiff;
    means fusion_id trt fusion_id*trt;
    /* Contrasts---------------------------------1---2---3---4---5---6---7---8---9---10---11---12-*/
    contrast 'AHCSxAHBerM' trt                   0   1  -1   0   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHBerF' trt             1   0   0  -1   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHCS' trt               0   0   1  -1   0   0   0   0   0    0    0    0 ;
    contrast 'AHBerMxAHBerF' trt                 1  -1   0   0   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSxAHBerF' trt                   1   0  -1   0   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHBerM' trt             0   1   0  -1   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSxAHFruP14440' trt              0   0   1   0   0   0   0  -1   0    0    0    0 ;
    contrast 'AHBerMxAHFruP14440' trt            0   1   0   0   0   0   0  -1   0    0    0    0 ;
    contrast 'AHCSxAHFruW12ChaM5' trt            0   0   1   0   0   0   0   0  -1    0    0    0 ;
    contrast 'AHBerMxAHFruW12ChaM5' trt          0   1   0   0   0   0   0   0  -1    0    0    0 ;
    contrast 'AHCSxAHMaleFruM(A)' trt            0   0   1   0   0   0   0   0   0   -1    0    0 ;
    contrast 'AHBerMxAHMaleFruM(A)' trt          0   1   0   0   0   0   0   0   0   -1    0    0 ;
    contrast 'AHCSxAHMaleFruM(B)' trt            0   0   1   0   0   0   0   0   0    0   -1    0 ;
    contrast 'AHBerMxAHMaleFruM(B)' trt          0   1   0   0   0   0   0   0   0    0   -1    0 ;
    contrast 'AHCSxAHMaleFruM(C)' trt            0   0   1   0   0   0   0   0   0    0    0   -1 ;
    contrast 'AHBerMxAHMaleFruM(C)' trt          0   1   0   0   0   0   0   0   0    0    0   -1 ;
    contrast 'AHCSFemalexAHFemaleFruM(A)' trt    0   0   0   1  -1   0   0   0   0    0    0    0 ;
    contrast 'AHBerFxAHFemaleFruM(A)' trt        1   0   0   0  -1   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
    contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHFemaleFruM(C)' trt    0   0   0   1   0   0  -1   0   0    0    0    0 ;
    contrast 'AHBerFxAHFemaleFruM(C)' trt        1   0   0   0   0   0  -1   0   0    0    0    0 ;

    output out=fru.resid_by_symbol r=resid p=pred;
    ods output modelanova=model_fru
               contrasts=fru.contrast_by_symbol
               lsmeans=fru.lsmeans_by_symbol
               means=fru.means_by_symbol;
    run;
    quit;

ods listing;

data fru.model_fru_by_symbol;
    set model_fru;
    where hypothesisType = 1;
    drop dependent hypothesisType ;
    run;

/* Check means to make sure they are the same
proc sort data=fru.all_coverage_counts_with_key;
    by trt fusion_id;
    run; * 5245230 = 87*60290;

data mean_check;
    set fru.all_coverage_counts_with_key;
    where fusion_id = "F10023_SI";
    run;

proc means data = mean_check noprint;
    class trt fusion_id;
    var logrpkm;
    output out=mean_check2 mean=mean;
    run;
*/

ods listing close;
proc glm data=fru.foranalysis_no_multi plots=none;
    by symbol_cat;
    class trt fusion_id;
    model logrpkm = fusion_id trt fusion_id*trt ;
    /* Contrasts---------------------------------1---2---3---4---5---6---7---8---9---10---11---12-*/
    contrast 'AHCSxAHBerM' trt                   0   1  -1   0   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHBerF' trt             1   0   0  -1   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHCS' trt               0   0   1  -1   0   0   0   0   0    0    0    0 ;
    contrast 'AHBerMxAHBerF' trt                 1  -1   0   0   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSxAHBerF' trt                   1   0  -1   0   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHBerM' trt             0   1   0  -1   0   0   0   0   0    0    0    0 ;
    contrast 'AHCSxAHFruP14440' trt              0   0   1   0   0   0   0  -1   0    0    0    0 ;
    contrast 'AHBerMxAHFruP14440' trt            0   1   0   0   0   0   0  -1   0    0    0    0 ;
    contrast 'AHCSxAHFruW12ChaM5' trt            0   0   1   0   0   0   0   0  -1    0    0    0 ;
    contrast 'AHBerMxAHFruW12ChaM5' trt          0   1   0   0   0   0   0   0  -1    0    0    0 ;
    contrast 'AHCSxAHMaleFruM(A)' trt            0   0   1   0   0   0   0   0   0   -1    0    0 ;
    contrast 'AHBerMxAHMaleFruM(A)' trt          0   1   0   0   0   0   0   0   0   -1    0    0 ;
    contrast 'AHCSxAHMaleFruM(B)' trt            0   0   1   0   0   0   0   0   0    0   -1    0 ;
    contrast 'AHBerMxAHMaleFruM(B)' trt          0   1   0   0   0   0   0   0   0    0   -1    0 ;
    contrast 'AHCSxAHMaleFruM(C)' trt            0   0   1   0   0   0   0   0   0    0    0   -1 ;
    contrast 'AHBerMxAHMaleFruM(C)' trt          0   1   0   0   0   0   0   0   0    0    0   -1 ;
    contrast 'AHCSFemalexAHFemaleFruM(A)' trt    0   0   0   1  -1   0   0   0   0    0    0    0 ;
    contrast 'AHBerFxAHFemaleFruM(A)' trt        1   0   0   0  -1   0   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
    contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
    contrast 'AHCSFemalexAHFemaleFruM(C)' trt    0   0   0   1   0   0  -1   0   0    0    0    0 ;
    contrast 'AHBerFxAHFemaleFruM(C)' trt        1   0   0   0   0   0  -1   0   0    0    0    0 ;

    output out=fru.resid_by_symbol r=resid p=pred;
    ods output modelanova=model_fru
               contrasts=fru.contrast_by_symbol
               ;
    run;
    quit;
ods listing;
