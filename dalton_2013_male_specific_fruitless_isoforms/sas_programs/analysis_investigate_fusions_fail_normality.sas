/********************************************************************************
* There were around 50 genes called out specifically in the manuscript. I
* looked at the results table and only saw 5 that were identified as being
* {ind,rep} in a treatment group because of a fusion that was flagged as
* failing normality.
*
* I am going to look more specifically at these fusions here.
********************************************************************************/

data subset;
    set FRU.subset_fail_normality;
    run;

    proc univariate data=S37849_SI  normal plot;
        var resid;
        run;
/* 
    AH_BerF
    AH_BerM
    AH_CS
    AH_CSFemale
    AH_Female_FruM(A
    AH_Female_FruM(B
    AH_Female_FruM(C
    AH_FruP14_440
    AH_FruW12_ChaM5
    AH_Male_FruM(A)
    AH_Male_FruM(B)
    AH_Male_FruM(C)

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
*/

/* S52794_SI: fru FruMB ind female */ 
* OK;

    proc glm data=subset; * Check and make sure contrast is significant so we know we picked the right one;
        where fusion_id="S52794_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
        output out=S52794_SI r=resid p=pred;
        run;

    proc gplot data=S52794_SI;
		by fusion_id;
		plot resid*pred=trt;
		run; * slight hetero;

    proc sort data=S52794_SI;
		by trt;
        run;

    proc univariate data=S52794_SI  normal plot;
		where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(B)" ;
		by trt;
		var resid;
		run; * AH_BerF shows a lot of variance in resid;
		
    proc univariate data=S52794_SI  normal plot;
		where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(B)" ;
		by trt;
		var logrpkm;
		run; * Clear Separation of trt from controls OK!;
		
    proc mixed data=subset ; 
        * Run model using proc mixed accounting for the differences in trt
        * variation. You expect the p-values to get more significant;
        where fusion_id="S52794_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
        repeated /group=trt;
        run; * OK;

/* S52798_SI: fru FruMB ind female */
* OK;

    proc glm data=subset; 
        * Check and make sure contrast is significant so we know we picked the
        * right one;
        where fusion_id="S52798_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
        output out=S52798_SI r=resid p=pred;
        run;

    proc gplot data=S52798_SI;
		by fusion_id;
		plot resid*pred=trt;
		run;* heterosch;

    proc sort data=S52798_SI;
		by trt;
        run;

    proc univariate data=S52798_SI  normal plot;
		where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(B)" ;
		by trt;
		var resid;
		run; * AH_BerF has a little more variance than others;
		
    proc univariate data=S52798_SI  normal plot;
		where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(B)" ;
		by trt;
		var logrpkm;
		run; * Clear separation, OK!;
		
    proc mixed data=subset ;
        * Run model using proc mixed accounting for the differences in trt
        * variation. You expect the p-values to get more significant;
        where fusion_id="S52798_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
        repeated /group=trt;
        run; * OK;

/* S47249_SI: dsx FruMC ind female */
* OK;

    proc glm data=subset ;
        * Check and make sure contrast is significant so we know we picked the
        * right one;
        where fusion_id="S47249_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(C)' trt    0   0   0   1   0   0  -1   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(C)' trt        1   0   0   0   0   0  -1   0   0    0    0    0 ;
        output out=S47249_SI r=resid p=pred;
        run;

    proc gplot data=S47249_SI;
        by fusion_id;
        plot resid*pred=trt;
        run; * bowtie shape resids;

    proc sort data=S47249_SI;
        by trt;
        run;

    proc univariate data=S47249_SI  normal plot;
        where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(C)" ;
        by trt;
        var resid;
        run; * AH_BerF has a lot more variance than others;

    proc univariate data=S47249_SI  normal plot;
        where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(C)" ;
        by trt;
        var logrpkm;
        run; * Separation, but close OK!;

    proc mixed data=subset ;
        * Run model using proc mixed accounting for the differences in trt
        * variation. You expect the p-values to get more significant;
        where fusion_id="S47249_SI "; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(C)' trt    0   0   0   1   0   0  -1   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(C)' trt        1   0   0   0   0   0  -1   0   0    0    0    0 ;
        repeated /group=trt;
        run; * Ok still significant;

/* S7129_SI: TOR FruMB ind male */
* WARNING;

    proc glm data=subset ;
        where fusion_id="S7129_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSxAHMaleFruM(B)' trt            0   0   1   0   0   0   0   0   0    0   -1    0 ;
        contrast 'AHBerMxAHMaleFruM(B)' trt          0   1   0   0   0   0   0   0   0    0   -1    0 ;
        output out=S7129_SI r=resid p=pred;
        run;

    proc gplot data=S7129_SI;
        by fusion_id;
        plot resid*pred=trt;
        run; * Major hetorsche;

    proc sort data=S7129_SI;
        by trt;
        run;

    proc univariate data=S7129_SI  normal plot;
        where trt="AH_CS" or trt="AH_BerM" or trt="AH_Male_FruM(B)" ;
        by trt;
        var resid;
        run; * AH_BerM has a lot more variance;

    proc univariate data=S7129_SI  normal plot;
        where trt="AH_CS" or trt="AH_BerM" or trt="AH_Male_FruM(B)" ;
        by trt;
        var logrpkm;
        run; * OK separation;

    proc mixed data=subset ;
        where fusion_id="S7129_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSxAHMaleFruM(B)' trt            0   0   1   0   0   0   0   0   0    0   -1    0 ;
        contrast 'AHBerMxAHMaleFruM(B)' trt          0   1   0   0   0   0   0   0   0    0   -1    0 ;
        repeated /group=trt;
        run; * Comparison with AH_BerM loses significance. WARNING!;

/* S37849_SI: Dh31-R1 FruMB ind male */
* WARNING!;

    proc glm data=subset ;
        where fusion_id="S37849_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSxAHMaleFruM(B)' trt            0   0   1   0   0   0   0   0   0    0   -1    0 ;
        contrast 'AHBerMxAHMaleFruM(B)' trt          0   1   0   0   0   0   0   0   0    0   -1    0 ;
        output out=S37849_SI r=resid p=pred;
        run; * P-value NS, I checked the FDR and it just passes our FDR criterion. It is ok but not great;

    proc gplot data=S37849_SI;
        by fusion_id;
        plot resid*pred=trt;
        run; * Lots of herterosche;

    proc sort data=S37849_SI;
        by trt;
        run;

    proc univariate data=S37849_SI  normal plot;
        where trt="AH_CS" or trt="AH_BerM" or trt="AH_Male_FruM(B)";
        by trt;
        var resid;
        run; * AH_CS has a lot of variance. Also AH_BerM looks to have little data;

    proc univariate data=S37849_SI normal plot;
        where trt="AH_CS" or trt="AH_BerM" or trt="AH_Male_FruM(B)";
        by trt;
        var logrpkm;
        run; * Dose not look that different WARNING!;

    proc mixed data=S37849_SI;
        by fusion_id;
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSxAHMaleFruM(B)' trt            0   0   1   0   0   0   0   0   0    0   -1    0 ;
        contrast 'AHBerMxAHMaleFruM(B)' trt          0   1   0   0   0   0   0   0   0    0   -1    0 ;
        repeated /group=trt;        
        run; * AH_CS comparison still is not significant WARNING!;
 
/* S38862_SI: Ir51a FruMA ind male, FruMC ind male */
* WARNING!;

    proc glm data=subset ;
        where fusion_id="S38862_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSxAHMaleFruM(A)' trt            0   0   1   0   0   0   0   0   0   -1    0    0 ;
        contrast 'AHBerMxAHMaleFruM(A)' trt          0   1   0   0   0   0   0   0   0   -1    0    0 ;
        contrast 'AHCSxAHMaleFruM(C)' trt            0   0   1   0   0   0   0   0   0    0    0   -1 ;
        contrast 'AHBerMxAHMaleFruM(C)' trt          0   1   0   0   0   0   0   0   0    0    0   -1 ;
        output out=S38862_SI r=resid p=pred;
        run; * AH_BerM vs FruMA is NS WARNING!;

    proc gplot data=S38862_SI;
        by fusion_id;
        plot resid*pred=trt;
        run; * Not too bad;

    proc sort data=S38862_SI;
        by trt;
        run;

    proc univariate data=S38862_SI  normal plot;
        where trt="AH_CS" or trt="AH_BerM" or trt="AH_Male_FruM(C)" or trt="AH_Male_FruM(A)" ;
        by trt;
        var resid;
        run; * AH_Male_FruM(C) has a lot of variation;

    proc univariate data=S38862_SI normal plot;
        where trt="AH_CS" or trt="AH_BerM" or trt="AH_Male_FruM(C)" or trt="AH_Male_FruM(A)" ;
        by trt;
        var logrpkm;
        run; * Clear separation OK;

    proc mixed data=S38862_SI;
        by fusion_id;
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSxAHMaleFruM(A)' trt            0   0   1   0   0   0   0   0   0   -1    0    0 ;
        contrast 'AHBerMxAHMaleFruM(A)' trt          0   1   0   0   0   0   0   0   0   -1    0    0 ;
        contrast 'AHCSxAHMaleFruM(C)' trt            0   0   1   0   0   0   0   0   0    0    0   -1 ;
        contrast 'AHBerMxAHMaleFruM(C)' trt          0   1   0   0   0   0   0   0   0    0    0   -1 ;
        repeated /group=trt;        
        run; * AH_FruMC is NS WARNING!;

/* S24566_SI: CG42330 (dscam4) FruMA ind female, FruMB ind female */
* OK;

    proc glm data=subset ;
        where fusion_id="S24566_SI"; 
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(A)' trt    0   0   0   1  -1   0   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(A)' trt        1   0   0   0  -1   0   0   0   0    0    0    0 ;
        contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
        output out=S24566_SI r=resid p=pred;
        run; * OK;

    proc gplot data=S24566_SI;
        by fusion_id;
        plot resid*pred=trt;
        run; * Not too bad;

    proc sort data=S24566_SI;
        by trt;
        run;

    proc univariate data=S24566_SI  normal plot;
        where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(A)" or trt="AH_Female_FruM(B)" ;
        by trt;
        var resid;
        run; * CS, A, and B have a lot of variation;

    proc univariate data=S24566_SI normal plot;
        where trt="AH_CSFemale" or trt="AH_BerF" or trt="AH_Female_FruM(A)" or trt="AH_Female_FruM(B)" ;
        by trt;
        var logrpkm;
        run; * Clear separation OK;

    proc mixed data=S24566_SI;
        by fusion_id;
        class trt ;
        model logrpkm = trt ;
        contrast 'AHCSFemalexAHFemaleFruM(A)' trt    0   0   0   1  -1   0   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(A)' trt        1   0   0   0  -1   0   0   0   0    0    0    0 ;
        contrast 'AHCSFemalexAHFemaleFruM(B)' trt    0   0   0   1   0  -1   0   0   0    0    0    0 ;
        contrast 'AHBerFxAHFemaleFruM(B)' trt        1   0   0   0   0  -1   0   0   0    0    0    0 ;
        repeated /group=trt;        
        run; * OK!;

