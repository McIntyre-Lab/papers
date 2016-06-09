/*******************************************************************************
* Filename: ase_summarize_ase_counts.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: Does a bunch of counts and summaries to get a better grasp on
* the ASE results.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Number of Genotypes Remaining */
    * 49 Genotypes remain after all filters;
    data geno;
        set CEGS.clean_ase_sbs;
        keep line;
        run;

    proc sort data=geno nodupkey;
        by line;
        run;

/* Number of exonic regions without AI */
    * 1301 exonic regions with no AI;
    * 4090 exonic regions with some AI;
    data noAI;
        set CEGS.clean_ase_sbs;
        margin_AI =  flag_AI_combined_m + flag_AI_combined_v;
        run;

    proc sort data=noAI;
        by fusion_id;
        run;

    proc means data=noAI noprint;
        by fusion_id;
        output out=sum sum(margin_AI)=sum_margin_AI;
        run;

    data flag_no_ai;
        set sum;
        if sum_margin_AI eq 0 then flag_no_ai = 1;
        else flag_no_ai = 0;
        keep fusion_id flag_no_ai;
        run;

    proc freq data=flag_no_ai;
        tables flag_no_ai;
        run;

    * Clean up;
    proc datasets nolist;
        delete sum;
        run; quit;

/* Drop exonic regions with no AI */
    proc sort data=noAI;
        by fusion_id;
        run;

    proc sort data=flag_no_ai;
        by fusion_id;
        run;

    data exonWithAI;
        merge noAI (in=in1) flag_no_ai (in=in2);
        by fusion_id;
        if in1 and flag_no_ai eq 0;
        if margin_AI > 0 then flag_AI_margin = 1;
        else flag_AI_margin = 0;
        drop flag_no_ai;
        run;

/* Number of genotypes with AI per exonic regions */
    * 14 exons have AI in more than 40 genotypes;
    * 31 exons show AI in 100% of genotypes tested;
    * 3004 exons show AI in 5 or fewer genotypes;
    * 570 exons show AI in 10% or fewer genotypes tested;
    proc freq data=exonWithAI noprint;
        by fusion_id;
        tables flag_AI_margin /out=freqs;
        run;

    data genoAiPerFus;
        set freqs;
        where flag_AI_margin eq 1;
        label count = ' ';
        label percent = ' ';
        keep fusion_id count percent;
        run;

    proc sort data=genoAiPerFus;
        by DESCENDING count;
        run;

    proc sort data=genoAiPerFus;
        by DESCENDING percent;
        run;

    proc sort data=genoAiPerFus;
        by count;
        run;

    proc sort data=genoAiPerFus;
        by percent;
        run;

    * Clean up;
    proc datasets nolist;
        delete freqs
        run; quit;

/* Number of exonic regions with AI per genotype */
    * NOTE: I have removed the 1301 exons that had no AI in any line. These
    * number are also based on mated and virgin combined.;
    * w47 had the most exons with AI (n=1244, ~50% of exons);
    * r21 had the fewest exons with AI (n=143, ~24% of exons);
    * r857 had the highest percent of exons with AI (~56% exons tested);
    * r380 had the lowest percent of exons with AI (~13% exons tested);
    * Most genotypes appear to have ~25-30% of remaining exons with AI;

    proc sort data=exonWithAI;
        by line;
        run;

    proc freq data=exonWithAI noprint;
        by line;
        tables flag_AI_margin /out=freqs;
        run;

    data exonAiPerGeno;
        set freqs;
        where flag_AI_margin eq 1;
        label count = ' ';
        label percent = ' ';
        keep line count percent;
        run;

    proc sort data=exonAiPerGeno;
        by DESCENDING count;
        run;

    proc sort data=exonAiPerGeno;
        by DESCENDING percent;
        run;

    proc sort data=exonAiPerGeno;
        by count;
        run;

    proc sort data=exonAiPerGeno;
        by percent;
        run;

    * Clean up;
    proc datasets nolist;
        delete freqs
        run; quit;

/* Mated and Virgin Genotype Summaries */
    * Most of the time, mated and virgin seem pretty similar. There are 10
    * genotypes that have more than a 100 exon difference between
    * mated and virgin.
    *       r502 r101 w68 r799 w64 r315 r340 w38 r427 w52
    ;

    proc sort data=exonWithAI;
        by line;
        run;

    proc means data=exonWithAI noprint;
        by line;
        output out=sums sum(flag_AI_combined_m)=mated_AI sum(flag_AI_combined_v)=virgin_AI;
        run;

    data mvGenoWithAI;
        set sums;
        keep line mated_AI virgin_AI;
        run;

    proc sort data=mvGenoWithAI;
        by mated_AI virgin_AI;
        run;

    data mv_big_diff;
        set mvGenoWithAI;
        if abs(mated_AI - virgin_AI) > 100;
        run;

    * Clean up;
    proc datasets nolist;
        delete sums;
        delete mv_big_diff;
        run; quit;

/* Mated and Virgin Genotype Summaries */
    * There are 469 exonic regions where mated showed no AI and virgin had at
    * least 1 genotype with AI.
    *
    * There are 491 exonic regions where virgin showed no AI and mated had at
    * least 1 genotype with AI.
    *
    * There is only 1 exonic region where mated and virgin are more than 10
    * genotypes different. [F9836_SI, M=14, V=27]
    ;
    proc sort data=exonWithAI;
        by fusion_id;
        run;

    proc means data=exonWithAI noprint;
        by fusion_id;
        output out=sums sum(flag_AI_combined_m)=mated_AI sum(flag_AI_combined_v)=virgin_AI;
        run;

    data mvFusWithAI;
        set sums;
        keep fusion_id mated_AI virgin_AI;
        run;

    proc sort data=mvFusWithAI;
        by mated_AI virgin_AI;
        run;

    data mv_big_diff;
        set mvFusWithAI;
        if abs(mated_AI - virgin_AI) > 10;
        run;

    data onOff;
        set mvFusWithAI;
        where (mated_AI eq 0 and virgin_AI gt 0) or (mated_AI gt 0 and virgin_AI eq 0);
        run;

    * Clean up;
    proc datasets nolist;
        delete sums;
        delete mv_big_diff;
        run; quit;

/* How many exons show extreme AI by genotype */
    data ind_extreme;
        set exonWithAI;
        if flag_AI_combined_m eq 1 then do;
            if q5_mean_theta_m < 0.4 then ind_mated_extreme = -1;
            else if q5_mean_theta_m > 0.6 then ind_mated_extreme = 1;
            else ind_mated_extreme = 0;
        end;
        else ind_mated_extreme = 0;
        if flag_AI_combined_v eq 1 then do;
            if q5_mean_theta_v < 0.4 then ind_virgin_extreme = -1;
            else if q5_mean_theta_v > 0.6 then ind_virgin_extreme = 1;
            else ind_virgin_extreme = 0;
        end;
        else ind_virgin_extreme = 0;
        keep line fusion_id ind_mated_extreme ind_virgin_extreme;
        run;

    proc sort data=ind_extreme;
        by line;
        run;

    proc freq data=ind_extreme noprint;
        by line;
        tables ind_mated_extreme /out=mfreq;
        tables ind_virgin_extreme /out=vfreq;
        run;

    proc transpose data=mfreq out=mflip;
        by line;
        var count;
        id ind_mated_extreme;
        run;

    data mated_geno_extreme_AI_table;
        set mflip;
        rename N1 = extreme_AI_toward_line;
        rename _0 = no_extreme_AI;
        rename _1 = extreme_AI_toward_tester;
        if N1 = . then N1 = 0;
        if _0 = . then _0 = 0;
        if _1 = . then _1 = 0;
        drop _name_ _label_;
        run;

    proc transpose data=vfreq out=vflip;
        by line;
        var count;
        id ind_virgin_extreme;
        run;

    data virgin_geno_extreme_AI_table;
        set vflip;
        rename N1 = extreme_AI_toward_line;
        rename _0 = no_extreme_AI;
        rename _1 = extreme_AI_toward_tester;
        if N1 = . then N1 = 0;
        if _0 = . then _0 = 0;
        if _1 = . then _1 = 0;
        drop _name_ _label_;
        run;

    * Clean up;
    proc datasets nolist;
        delete mfreq;
        delete vfreq;
        delete mflip;
        delete vflip;
        run; quit;

/* How many exons show extreme AI by exon */
    proc sort data=ind_extreme;
        by fusion_id;
        run;

    proc freq data=ind_extreme noprint;
        by fusion_id;
        tables ind_mated_extreme /out=mfreq;
        tables ind_virgin_extreme /out=vfreq;
        run;

    proc transpose data=mfreq out=mflip;
        by fusion_id;
        var count;
        id ind_mated_extreme;
        run;

    data mated_fus_extreme_AI_table;
        set mflip;
        rename N1 = extreme_AI_toward_line;
        rename _0 = no_extreme_AI;
        rename _1 = extreme_AI_toward_tester;
        if N1 = . then N1 = 0;
        if _0 = . then _0 = 0;
        if _1 = . then _1 = 0;
        drop _name_ _label_;
        run;

    proc transpose data=vfreq out=vflip;
        by fusion_id;
        var count;
        id ind_virgin_extreme;
        run;

    data virgin_fus_exteme_AI_table;
        set vflip;
        rename N1 = extreme_AI_toward_line;
        rename _0 = no_extreme_AI;
        rename _1 = extreme_AI_toward_tester;
        if N1 = . then N1 = 0;
        if _0 = . then _0 = 0;
        if _1 = . then _1 = 0;
        drop _name_ _label_;
        run;

    * Clean up;
    proc datasets nolist;
        delete mfreq;
        delete vfreq;
        delete mflip;
        delete vflip;
        run; quit;
