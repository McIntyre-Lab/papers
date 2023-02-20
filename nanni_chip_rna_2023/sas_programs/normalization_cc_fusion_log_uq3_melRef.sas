
libname ortho "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/sas_data/RNAseq_ortho";
libname drosRNA "!MCLAB/Dros_CHIP_RNA_ms/sas_data/RNAseq";
libname df  "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/svn_lmm_dros_head_data/ethanol_srna/sas_data";


filename mymacros "!MCLAB/maize_ozone_final/2014/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);




/********************************************************************************
* Normalization using the log2 of the upper quartile adjusted APN.
* into normalization - fusions where flag_sim_f_on = 1 AND flag_mel_f_on = 1
********************************************************************************/

/* Create clean dataset 

drop sim_12_m_noEtoh_rep2 before norm

UQmult of 57 = q3 across all samples of the median of the reads in region

Also calculating fudge factor based on mapped reads

        mapped_FF = (sum_mapped / (median*10**5))
        map_apn = (apn*mapped_FF); 
        log_map_apn = log(map_apn);
        uq_apn = (apn/q3)*57;       
        uq_ff = 57/q3;
        log_uq_apn = log(UQ_apn);

*/



data dsgn_file ;
retain sampleID ;
length sampleID $21.;
set df.design_ethanol_srna_rnaseq ;
sampleID = compress(species||'_'||genotype||'_'||sex||'_'||treatment||'_rep'||arbitrary_rep);
rep = compress('rep'||arbitrary_rep);
keep sampleID species genotype sex rep treatment;
run ;

proc sort data = dsgn_file ;
by sampleID ;
proc sort data=ortho.cvrg_fusion_melref_stack;
by sampleID ;
run ;

data clean1;
merge ortho.cvrg_fusion_melref_stack(in=in1) dsgn_file (in=in2);
by sampleID ;
run;

proc sort data = clean1 ;
by featureID ;
proc sort data=ortho.oncall_sample_apn0_melRef;
by featureID;
run;

data clean2;
merge clean1 (in=in1) ortho.oncall_sample_apn0_melRef (in=in2);
by featureID ;
run;

data clean3 ;
set clean2 ;
if flag_F_on0 = 1 or flag_M_on0 = 1 ;
if sampleID = "sim_12_m_noEtoh_rep2"  then delete ;
run ;

/* Calculate With-in sample basic statistics */
    proc means data=clean3 noprint;
        class species genotype sex treatment rep;
        var reads_in_region;
        output out=mapped_reads sum=sum_mapped q3=q3 median=median;
        run;

    data totals;
        set mapped_reads;
        where _type_ = 31;
        keep species genotype sex treatment rep sum_mapped q3 median;
        run;

    /* Summarize statistics to experiment level */
        proc sort data=totals;
            by species sex ;
            run;

        proc means data=totals median ;
            var sum_mapped q3 median;
            run;
            /* q3 = 57.13  
            FF = (57 / sum_mapped)  */

        /* check q3 by sample -- NOTE: sim_12_m_noEtoh_rep2 q3 = 49, this is NOT low like it was in the RNAseq analysis, removing does NOT affect overall q3 value */
        /*proc sort data = totals; 
        by species genotype sex treatment rep ;
        run ;

        proc means data=totals median ;
            by species genotype sex treatment rep  ;
            var sum_mapped q3 median;
            output out = check_q3 mean = ;
            run;

            proc sort data = check_q3 ;
                by q3 ;
                run ;
     */
    /* Merge statistics onto dataset */
    proc sort data=clean3;
        by species genotype sex treatment rep;
        run;

    proc sort data=totals;
        by species genotype sex treatment rep;
        run;

    data stack oops;
        merge clean3 (in=in1) totals (in=in2);
        by species genotype sex treatment rep;
        if in1 and in2 then output stack;
        else output oops; * 0 obs yay!;
        run;


    /* NOTE MAKE SURE TO CHANGE THE UQ MULTIPLIER BY THE NEW Q3 

        mapped_FF = (sum_mapped / (median*10**5))
        map_apn = (apn*mapped_FF); 
        log_map_apn = log(map_apn);
        uq_apn = (apn/q3)*57;       
        uq_ff = 57/q3;
        log_uq_apn = log(UQ_apn);

*/
    data stack_uq;
        retain featureID;
        set stack;
        mapped_FF = (sum_mapped / (median*10**5)) ;
        map_apn = (apn*mapped_FF); 
        log_map_apn = log(map_apn); 
        uq_apn = (apn/q3)*57;
        uq_ff = 57/q3;
        log_uq_apn = log2(UQ_apn);
        run;

    proc freq data=stack_uq noprint;
        table featureID/ out=freqs;
        run;

    /* make perm */
data ortho.norm_stats_fusion_melRef ;
set stack_uq ;
run;


