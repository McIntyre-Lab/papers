/*******************************************************************************
* Filename: cis_eq_data_normalization.sas
*
* Author: Justin M Fear | jfear@ufl.edu
*
* Description: For the maren equations, I am also going to drop exonic regions
* with less than 10 genotypes. The maren equations make some assumptions about
* the population level sums. Obvisouly the more genotypes that are present for
* each fusions the better, but I am comfortable with as few as 10 genotypes.
*
*******************************************************************************/

/* Libraries
libname cegs '!MCLAB/cegs_ase_paper/sas_data';
libname DMEL551 '!MCLAB/useful_dmel_data/flybase551/sasdata';
libname GENELIST '!MCLAB/useful_dmel_data/gene_lists/sas_data';
filename mymacros '!MCLAB/cegs_ase_paper/sas_programs/mclib_SAS';
options SASAUTOS=(sasautos mymacros);
*/

/* Count the number of lines present for each fusion */
    proc sort data=CEGS.clean_ase_stack;
        by mating_status;
        run;

    proc freq data=CEGS.clean_ase_stack noprint;
        by mating_status;
        tables fusion_id /out=freqs;
        run;

    data cleanFreq;
        set freqs;
        if count ge 10;
        label count = ' ';
        rename count = geno_count;
        keep mating_status fusion_id count;
        run;

/* Merge and drop fusions with few genotypes */
    proc sort data=CEGS.clean_ase_stack;
        by mating_status fusion_id;
        run;

    proc sort data=cleanFreq;
        by mating_status fusion_id;
        run;

    data clean;
        merge  CEGS.clean_ase_stack (in=in1) cleanFreq (in=in2);
        by mating_status fusion_id;
        if in2;
        line_tester_avg = (sum_line + sum_tester) / 2;
        run; * 131700 obs;

/* Calculate Average Expression or mu in the maren equations */
    proc sort data=clean;
        by fusion_id mating_status line;
        run;

    proc means data=clean noprint;
        by fusion_id mating_status;
        output out=means mean(line_tester_avg)=mu;
        run;

/* Merge Mu to clean */
    proc sort data=clean;
        by fusion_id mating_status;
        run;

    proc sort data=means;
        by fusion_id mating_status;
        run;

    data muclean;
        merge  clean (in=in1) means (in=in2 drop=_type_ _freq_);
        by fusion_id mating_status;
        run;

/* Coverage Center Data */
    data CEGS.data_for_cis_eq;
        set muclean;
        line_prop_adj = sum_line / (sum_line + sum_tester) * 1000;
        tester_prop_adj = sum_tester / (sum_line + sum_tester) * 1000;
        line_mean_center = sum_line - mu;
        tester_mean_center = sum_tester - mu;
        keep line mating_status fusion_id q5_mean_theta flag_AI_combined
        mean_apn line_prop_adj tester_prop_adj line_mean_center tester_mean_center geno_count;
        run;

/* Clean Up */
    proc datasets nolist;
        delete clean;
        delete cleanFreq;
        delete freqs;
        delete means;
        delete muclean;
        run; quit;
