/*******************************************************************************
* Filename: ase_summarize_ase_filters.sas
*
* Author: Justin M Fear | jfear@ufl.edu 
*
* Description: I have explored the ASE results and have found some concerns.
* After a lot of trail and error I have come up with a set filters to apply
* before proceeding with the analysis. 
*
* For a list of all of the things I have tried see:
* ../scripts/00_Notebook_TOC.ipynb
*
* Here I implement the filtering strategy that I describe here:
* ../scripts/ase_summary/ase_filters.ipynb
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Little macro to print how many fusions there are */
    %macro fusionCnt(dat);
        data tmp;
            set &dat;
            keep fusion_id;
            run;

        proc sort data=tmp nodupkey;
            by fusion_id;
            run;

        data _null_;
            put nobs=;
            stop;
            set tmp nobs=nobs;
            run;

        * Clean up;
        proc datasets nolist;
            delete tmp;
            run; quit;
    %mend;

/* Drop exonic regions flagged in 100 genome simulation */
    proc sort data=CEGS.qsim_emp_theta_w_flag;
        by fusion_id;
        run;

    %fusionCnt(CEGS.qsim_emp_theta_w_flag); * 49,947 obs;

    proc sort data=CEGS.exon_drop_list_100_genome;
        by fusion_id;
        run;

    data clean100Genome;
        merge CEGS.qsim_emp_theta_w_flag (in=in1) CEGS.exon_drop_list_100_genome (in=in2);
        by fusion_id;
        if in1 and not in2;
        run;

    %fusionCnt(clean100Genome); * 49,817 obs;

/* Drop exonic regions that have an APN < 25 */
    data goodAPN;
        set clean100Genome;
        if mean_apn < 25 then delete;
        run;

    %fusionCnt(goodAPN); * 18,026 obs;

/* Flag and drop exonic regions not in at least 10% genotypes */
    proc freq data=goodAPN noprint;
        tables line /out=lineFreqs;
        run; * 68 genotypes;

    proc sort data=goodAPN;
        by mating_status;
        run;

    proc freq data=goodAPN noprint;
        by mating_status;
        tables fusion_id /out=freqs;
        run;

    data flag_low_per_geno;
        set freqs;
        if count < 68 * .10 then flag_low_per_geno = 1;
        else flag_low_per_geno = 0;
        keep mating_status fusion_id flag_low_per_geno;
        run;

    proc datasets nolist;
        delete freqs lineFreqs;
        run; quit;

    proc sort data=goodAPN;
        by mating_status fusion_id;
        run;

    data goodPerGeno;
        merge goodAPN (in=in1) flag_low_per_geno (in=in2);
        by mating_status fusion_id;
        if in1 and flag_low_per_geno eq 0;
        drop flag_low_per_geno;
        run;

    %fusionCnt(goodPerGeno); * 6,644 obs;

/* Flag regions not in both environments */
    proc sort data=goodPerGeno;
        by line fusion_id;
        run;

    proc freq data=goodPerGeno noprint;
        by line ;
        tables fusion_id /out=freqs;
        run;

    data flag_one_environment;
        set freqs;
        if count eq 1 then flag_one_environment = 1;
        else if count eq 2 then flag_one_environment = 0;
        else flag_one_environment = -999;
        keep line fusion_id flag_one_environment;
        run;

    proc datasets nolist;
        delete freqs;
        run; quit;

    proc sort data=goodPerGeno;
        by line fusion_id;
        run;

    data goodEnv;
        merge goodPerGeno (in=in1) flag_one_environment (in=in2);
        by line fusion_id;
        if in1 and flag_one_environment eq 0;
        drop flag_one_environment;
        run;

    %fusionCnt(goodEnv); * 5,391 obs;

/* Flag and drop genotypes that show bias */
    * For a given line*mating_status we expect that most of the exons will not
    * show bias. This would result in a point mass around 0.5. I want to drop
    * genotypes whose point mass (median) is less than 0.4 or greater than 0.6.
    ;
    proc sort data=goodEnv;
        by line mating_status;
        run;

    proc means data=goodEnv;
        by line mating_status;
        output out=median median(q5_mean_theta)=med;
        run;

    data flag_geno_bias;
        set median;
        if med <= 0.4 or med >= 0.6 then flag_geno_bias = 1;
        else flag_geno_bias = 0;
        run;

    proc means data=flag_geno_bias noprint;
        by line;
        output out=sum sum(flag_geno_bias)=ind_geno_bias;
        run;

    data ind_geno_bias;
        set sum;
        keep line ind_geno_bias;
        run;

    data goodGenoBias;
        merge goodEnv (in=in1) ind_geno_bias (in=in2);
        by line;
        if in1 and ind_geno_bias eq 0;
        drop ind_geno_bias;
        run;

    %fusionCnt(goodGenoBias); * 5,391 obs;

    * Clean up;
    proc datasets nolist;
        delete median sum;
        run; quit;

/* Flag and drop genotypes with too few exonic regions */
    proc freq data=goodGenoBias noprint;
        by line;
        tables mating_status /out=freqs;
        run;

    data flag_low_exon_cnt;
        set freqs;
        if count <= 500 then flag_low_exon_cnt = 1;
        else flag_low_exon_cnt = 0;
        keep line mating_status flag_low_exon_cnt;
        run; 
        * Mated and virgin are the same, because I have already removed exons
        * that are not in both;

    data goodExonCnt;
        merge goodGenoBias (in=in1) flag_low_exon_cnt (in=in2);
        by line mating_status;
        if in1 and flag_low_exon_cnt = 0;
        drop flag_low_exon_cnt;
        run;

    %fusionCnt(goodExonCnt); * 5,391 obs;

    * Clean up;
    proc datasets nolist;
        delete freqs;
        run; quit;

/* Create permanant stack dataset */
    data CEGS.clean_ase_stack;
        set goodExonCnt;
        run;

/* Make M/V side-by-side */
    data mvstack;
        set goodExonCnt;
        keep line mating_status fusion_id q5_mean_theta flag_AI_combined;
        run;

    proc sort data=mvstack;
        by line fusion_id;
        run;

    data mvsbs;
        merge mvstack (in=in1 where=(mating_status = 'M') 
                       rename=(q5_mean_theta = q5_mean_theta_m
                               flag_AI_combined = flag_AI_combined_m)) 
              mvstack (in=in2 where=(mating_status = 'V') 
                       rename=(q5_mean_theta = q5_mean_theta_v
                       flag_AI_combined = flag_AI_combined_v));
        by line fusion_id;
        drop mating_status;
        run;

/* Create permanant sbs dataset */
    data CEGS.clean_ase_sbs;
        set mvsbs;
        run;

/* Merge on Gene information and export */
proc sort data=CEGS.clean_ase_sbs;
    by fusion_id;
    run;

proc sort data=DMEL551.fb551_si_fusions_unique;
    by fusion_id;
    run;

data merged;
    retain line fusion_id q5_mean_theta_m q5_mean_theta_v flag_AI_combined_m
    flag_AI_combined_v chrom start end symbol_cat fbgn_cat
    flag_AI_combined_m_and_v flag_AI_combined_m_or_v flag_m_and_v_opposite
    flag_discordant ;

    merge CEGS.clean_ase_sbs (in=in1)  DMEL551.fb551_si_fusions_unique (in=in2 keep=fusion_id symbol_cat fbgn_cat chrom start end);
    by fusion_id;
    if in1;

    if flag_AI_combined_m eq 1 and flag_AI_combined_v eq 1 then flag_AI_combined_m_and_v = 1;
    else flag_AI_combined_m_and_v = 0;

    if (flag_AI_combined_m eq 1 and flag_AI_combined_v eq 0) or (flag_AI_combined_m eq 0 and flag_AI_combined_v eq 1)then flag_AI_combined_m_or_v = 1;
    else flag_AI_combined_m_or_v = 0;

    if (q5_mean_theta_m > 0.5 and q5_mean_theta_v < 0.5) or (q5_mean_theta_m < 0.5 and q5_mean_theta_v > 0.5) then flag_m_and_v_opposite = 1;
    else flag_m_and_v_opposite = 0;

    if flag_AI_combined_m_and_v eq 1 and flag_m_and_v_opposite eq 1 then flag_discordant = 1;
    else flag_discordant = 0;

    run;

proc export data=merged outfile='!MCLAB/cegs_ase_paper/pipeline_output/ase_summary/exonic_regions_with_ai.csv' dbms=csv replace;
    putnames=yes;
    run;

/* Clean up */
    proc datasets nolist;
        delete clean100Genome;
        delete flag_low_per_geno;
        delete flag_one_environment;
        delete flag_geno_bias;
        delete ind_geno_bias;
        delete flag_low_exon_cnt;
        delete goodapn;
        delete goodenv;
        delete goodpergeno;
        delete goodGenoBias;
        delete goodExonCnt;
        delete mvsbs;
        delete mvstack;
        delete merged;
        run; quit;
