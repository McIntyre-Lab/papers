libname dros "!MCLAB/Dros_PB_ChIP/sasdata/RNAseq";


/*
* fusions

UQMult = 63!!

* for noEtoh and etoh:
    * Normalization using the log2 of the upper quartile adjusted APN. 
    * Only included fusions where flag_fusion_F_on = 1 and flag_fusion_M_on = 1 for APN gt 0.
    * Normalization is performed separately for Female and Male.

    input
        dros.df_pbrna_4_analysis
        dros.cnts_mel_fusions_apn_stack 
        dros.cnts_sim_fusions_apn_stack 
        dros.onCalls_fsn_&species._&trt._gt_apn0

    output
        dros.norm_stats_fsn_&species._F
        dros.norm_stats_fsn_&species._M

*/

%macro norm (species) ;


/* merge in design file */
proc sort data = dros.df_pbrna_4_analysis ;
by sampleID ;
proc sort data = dros.cnts_&species._fusions_apn_stack ;
by sampleID ;
run;

data &species._fusions ouch_&species;
merge dros.cnts_&species._fusions_apn_stack (in=in1) dros.df_pbrna_4_analysis (in=in2) ;
by sampleID ;
if in1 and in2 then output &species._fusions ;
else output ouch_&species. ;
run;


/* Create clean dataset - drop fusions with little expression (keep flag_fusion_F_on = 1 and flag_fusion_M_on = 1) */
proc sort data=&species._fusions;
by featureID ;
proc sort data = dros.onCalls_fsn_&species._etoh_gt_apn0 ;
by featureID ;

proc sort data = dros.onCalls_fsn_&species._noEtoh_gt_apn0 ;
by featureID ;
run;

data &species._clean1;
merge &species._fusions (in=in1) dros.onCalls_fsn_&species._etoh_gt_apn0  ;
by featureID ;
if flag_fusion_F_on = 1 and flag_fusion_M_on = 1;
drop flag_: ;
run;  

data &species._clean ;
merge &species._clean1 (in=in1) dros.onCalls_fsn_&species._noEtoh_gt_apn0 ;
by featureID ;
if flag_fusion_F_on = 1 and flag_fusion_M_on = 1;
drop flag_: ;
run;    


/* Calculate With-in sample basic statistics */
    proc means data = &species._clean noprint;
        class genotype sex rep;
        var reads_in_region;
        output out=&species._mapped_reads sum=sum_mapped q3=q3 median=median;
        run;

    data totals_&species.;
        set &species._mapped_reads;
        where _type_ = 7;
        keep genotype sex rep sum_mapped q3 median;
        run;

    /* Summarize statistics to experiment level */
        proc sort data=totals_&species.;
            by sex ;
            run;

        proc means data=totals_&species. median ;
            by sex;
            var sum_mapped q3 median;
            output out = exp_&species. q3=;
            run;

    /* Merge statistics onto dataset */
    proc sort data=&species._clean;
        by genotype sex rep;
    proc sort data=totals_&species. ;
        by genotype sex rep;
        run;

    data stack_&species. oops_&species.;
        merge &species._clean (in=in1) totals_&species. (in=in2);
        by genotype sex rep;
        if in1 and in2 then output stack_&species.;
        else output oops_&species.; * 0 obs yay!;
        run;
%mend ;

%norm (mel) ;
%norm (sim) ;


%macro uq_mult (species, sex) ;

    /* UQ MULTIPLIER = 63!!  */
    data stack_&species._&sex;
        retain featureID;
        set stack_&species. (where=(sex = "&sex"));
        uq_apn = (apn/q3)*63;       
        uq_ff = 63/q3;
        log_uq_apn = log2(UQ_apn);
        /* log_uq_apn = log2(UQ_apn + 1);   -- use UQ_apn + 1 when plotting and want to have 0's instead of missing values */
        run;


    /* Make sure that all fusions are present where flag_fusion_F_on = 1 and flag_fusion_M_on =1 */
    proc freq data=stack_&species._&sex noprint;
        table featureID/ out=freqs_&sex._&species.;

    /* make perm */
    data dros.norm_stats_fsn_&species._&sex.;
        set stack_&species._&sex;
        run;

%mend ;

/* uq_mult values come from q3 - mean cell in exp_&species._&trt datasets  */
%uq_mult (mel, f) ;
%uq_mult (mel, m) ;

%uq_mult (sim, f) ;
%uq_mult (sim, m) ;




