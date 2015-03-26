/*
 * REVISIONS: 12/21/2011
 *            - From Import I changed all_coverage_counts to
 *              fru.all_coverage_counts
 */

libname fru '!MCLAB/Fru_network/sasdata';

proc sort data=fru.all_coverage_counts;
    by sample_id;
    run;

proc sort data=fru.design_file;
    by sample_id;
    run;
        
data fru.all_coverage_counts_with_key;
    merge fru.all_coverage_counts (in=in1) fru.design_file (in=in2);
    by sample_id;
    logrpkm = log(rpkm);
    if in1;
    run; *3617460 obs;
    /* 
     * There are some going to be some from all_coverage_counts_with_key that
     * will not match anything in fru.design_file. This is because "Justin
     * Dalton" left samples out of the comparison for 07/05/2011 lanes 5,6,7,8.
     * Also 05/03/2011 (aka 04/26/2011) lane 8 was a bad sample.
    */

