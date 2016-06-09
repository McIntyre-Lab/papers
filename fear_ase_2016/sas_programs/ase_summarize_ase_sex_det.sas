/*******************************************************************************
* Filename: ase_summarize_ase_sex_det.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Look only at the sex determination genes and determine what genes
* have evidence for AI.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Merge on Sex Determiniation Genes */
    data WORK.SEX_DET;
        infile '!MCLAB/useful_dmel_data/gene_lists/pathways/sex_det.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
        informat symbol $43. ;
        informat primary_fbgn $11. ;
        format symbol $43. ;
        format primary_fbgn $11. ;
        input
            symbol $
            primary_fbgn $
            ;
        run;

    proc sort data=sex_det;
        by primary_fbgn;
        run;

    proc sort data=DMEL551.FB551_si_fusion_2_gene_id;
        by FBgn;
        run;

    data sex2fus oops;
        merge sex_det (in=in1) DMEL551.FB551_si_fusion_2_gene_id (in=in2 rename=(FBgn=primary_fbgn));
        by primary_fbgn;
        if in1 and not in2 then output oops;
        if in1 and in2 then output sex2fus;
        keep symbol primary_fbgn fusion_id;
        run; * oops had none!!;

    * Clean up;
    proc datasets nolist;
        delete oops sex_det;
        run; quit;

/* Merge Clean Data to Sex Determination Genes */
    proc sort data=sex2fus;
        by fusion_id;
        run;

    proc sort data=CEGS.clean_ase_sbs;
        by fusion_id;
        run;

    data merged;
        merge sex2fus (in=in1) CEGS.clean_ase_sbs (in=in2);
        by fusion_id;
        if in1;
        keep symbol primary_fbgn fusion_id line flag_AI_combined_m flag_AI_combined_v;
        run;

/* Freq of sex det fusions in Clean Data */
    * For each sex det gene, look at the number of fusions that are in the clean data.
    * Then figure out the frequency of times a gene is on in the clean data.
    ;
    proc sort data=merged;
        by symbol;
        run;

    proc freq data=merged noprint;
        by symbol;
        tables flag_AI_combined_m /out=mfreq missing;
        tables flag_AI_combined_v /out=vfreq missing;
        run;

    /* Mated */
        data mfreq2;
            set mfreq;
            if flag_AI_combined_m = 0 then AI_combined = 'per_on_no_AI_m' ;
            if flag_AI_combined_m = 1 then AI_combined = 'per_on_AI_m' ;
            if flag_AI_combined_m = . then AI_combined = 'per_off_m' ;
            keep line symbol percent AI_combined;
            run;

        proc transpose data=mfreq2 out=mflip;
            by symbol;
            var percent;
            id AI_combined;
            run;

        data mOff;
            set mflip;
            if per_off_m = . then per_off_m = 0;
            if per_on_no_AI_m = . then per_on_no_AI_m = 0;
            if per_on_AI_m = . then per_on_AI_m = 0;
            drop _name_ _label_;
            run;

    /* Virgin */
        data vfreq2;
            set vfreq;
            if flag_AI_combined_v = 0 then AI_combined = 'per_on_no_AI_v' ;
            if flag_AI_combined_v = 1 then AI_combined = 'per_on_AI_v' ;
            if flag_AI_combined_v = . then AI_combined = 'per_off_v' ;
            keep line symbol percent AI_combined;
            run;

        proc transpose data=vfreq2 out=vflip;
            by symbol;
            var percent;
            id AI_combined;
            run;

        data vOff;
            set vflip;
            if per_off_v = . then per_off_v = 0;
            if per_on_no_AI_v = . then per_on_no_AI_v = 0;
            if per_on_AI_v = . then per_on_AI_v = 0;
            drop _name_ _label_;
            run;

    /* Merge On Off Results */
        proc sort data=moff;
            by symbol;
            run;

        proc sort data=voff;
            by symbol;
            run;

        data off_on;
            merge moff (in=in1) voff (in=in2);
            by symbol;
            run;

        data CEGS.sex_det_on_off;
            set off_on;
            run;

    * Clean up;
    proc datasets nolist;
        delete mfreq vfreq;
        delete mfreq2 vfreq2;
        delete mflip vflip;
        delete mOff vOff;
        delete off_on;
        run; quit;

/* Summarize fusions to Gene Level */
    proc sort data=merged;
        by symbol line;
        run;

    proc means data=merged noprint;
        by symbol line;
        where line ne '' ;
        output out=sums sum(flag_AI_combined_m)= sum(flag_AI_combined_v)= /autoname;
        run;

    data sums2;
        set sums;
        if flag_AI_combined_m_Sum eq 0 then AI_mated = 'mated_geno_no_AI';
        else AI_mated = 'mated_geno_w_AI';
        if flag_AI_combined_v_Sum eq 0 then AI_virgin = 'virgin_geno_no_AI';
        else AI_virgin = 'virgin_geno_w_AI';
        keep symbol line AI_mated AI_virgin;
        run;
        
    proc freq data=sums2 noprint;
        by symbol;
        tables AI_mated / out=mfreq;
        tables AI_virgin / out=vfreq;
        run;

    proc transpose data=mfreq out=mflip;
        by symbol;
        var count;
        id AI_mated;
        run;

    proc transpose data=vfreq out=vflip;
        by symbol;
        var count;
        id AI_virgin;
        run;

    data mvAI;
        merge mflip (in=in1) vflip (in=in2);
        by symbol;
        if mated_geno_no_AI = . then mated_geno_no_AI = 0;
        if mated_geno_w_AI = . then mated_geno_w_AI = 0;
        if virgin_geno_no_AI = . then virgin_geno_no_AI = 0;
        if virgin_geno_w_AI = . then virgin_geno_w_AI = 0;
        drop _name_ _label_;
        run;

    data CEGS.sex_det_geno_AI;
        set mvAI;
        run;

    * Clean up;
    proc datasets nolist;
        delete sums sums2;
        delete mfreq vfreq;
        delete mflip vflip;
        run; quit;

/* Clean up */
    proc datasets nolist;
        delete merged;
        delete mvAI;
        delete sex2fus;
        run; quit;
