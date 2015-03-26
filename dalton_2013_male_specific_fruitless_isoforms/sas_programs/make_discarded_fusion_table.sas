libname fru '!MCLAB/arbeitman_fru_network/sasdata';
libname dmel530 '!MCLAB/useful_dmel_data/flybase530/sasdata';

data fru.flag_multigene;
    set dmel530.fb530_si_fusions_unique_flagged;
    if Genes_per_fusion > 1 then flag_multigene = 1; else flag_multigene = 0;
    keep fusion_id flag_multigene;
    run;

proc freq data=fru.flag_multigene;
    table flag_multigene;
    run;

/*
    flag_multigene    Frequency
    ---------------------------
                 0       57962
                 1        2329
*/

proc sort data=FRU.all_coverage_counts_with_key;
    by fusion_id;
    run;

proc transpose data=FRU.all_coverage_counts_with_key out=flip;
    var apn;
    by fusion_id;
    id sample_id;
    run;

data row_sums;
    set flip;
    row_sum = sum(_2011_07_05_5_ATCACG, _2011_07_05_5_CAGATC,
    _2011_07_05_5_TTAGGC, _2011_07_05_6_ACTTGA, _2011_07_05_6_CGATGT,
    _2011_07_05_6_TAGCTT, _2011_07_05_6_TGACCA, _2011_07_05_7_ACAGTG,
    _2011_07_05_7_GATCAG, _2011_07_05_8_ACAGTG, _2011_07_05_8_ATCACG,
    _2011_07_05_8_CGATGT, _2011_07_05_8_GCCAAT, _2011_05_03_2_GCCAAT,
    _2011_05_03_3_ACAGTG, _2011_05_03_5_TGACCA, _2011_05_03_1_ATCACG,
    _2011_05_03_5_TTAGGC, _2011_05_03_6_CGATGT, _2011_05_03_1_ACTTGA,
    _2011_05_03_2_ACTTGA, _2011_05_03_3_ACTTGA, _2011_05_03_5_ACTTGA,
    _2011_07_05_7_ACTTGA, _2011_07_05_7_CAGATC, _2011_07_05_8_GATCAG,
    _2011_07_05_1_TGACCA, _2011_07_05_2_ACAGTG, _2011_07_05_3_GCCAAT,
    _2011_07_05_1_GCCAAT, _2011_07_05_2_TGACCA, _2011_07_05_3_ACAGTG,
    _2011_07_05_5_CGATGT, _2011_07_05_5_TAGCTT, _2011_07_05_6_ATCACG,
    _2011_05_03_1_GATCAG, _2011_05_03_2_GATCAG, _2011_05_03_3_GATCAG,
    _2011_05_03_1_TAGCTT, _2011_05_03_2_TAGCTT, _2011_05_03_3_TAGCTT,
    _2011_05_03_7_TAGCTT, _2011_07_05_1_ATCACG, _2011_07_05_2_CGATGT,
    _2011_07_05_3_TTAGGC, _2011_07_05_1_TTAGGC, _2011_07_05_2_ATCACG,
    _2011_07_05_3_CGATGT, _2011_07_05_1_CAGATC, _2011_07_05_2_ACTTGA,
    _2011_07_05_3_GATCAG, _2011_05_03_1_CAGATC, _2011_05_03_5_GATCAG,
    _2011_05_03_6_ACTTGA, _2011_05_03_2_ATCACG, _2011_05_03_5_CGATGT,
    _2011_05_03_6_TAGCTT, _2011_05_03_3_TGACCA, _2011_05_03_6_ACAGTG,
    _2011_05_03_6_TTAGGC);
    run;

data fru.flag_no_reads;
    set row_sums;
    if row_sum = 0 then flag_no_reads = 1; else flag_no_reads = 0;
    keep fusion_id flag_no_reads;
    run;


proc sort data=fru.flag_multigene;
    by fusion_id;
    run;

proc sort data=fru.flag_no_reads;
    by fusion_id;
    run;

proc sort data=fru.flag_no_var;
    by fusion_id;
    run;

data merged;
    merge fru.flag_multigene fru.flag_no_reads fru.flag_no_var;
    by fusion_id;
    run;

proc freq data=merged;
    tables flag_multigene*flag_no_reads*flag_no_var;
    run;

proc freq data=merged(where=(flag_multigene = 0));
    table flag_no_reads;
    run;
    
/*
    flag_no_reads    Frequency
    --------------------------
                0       56443
                1        1519
*/

proc freq data=merged(where=(flag_multigene = 0 and flag_no_reads = 0));
    table flag_no_var;
    run;
    
/*
    flag_no_var    Frequency
    -------------------------
              0       45339
              1       11104
*/

/* create output table of these 11104 fusions */

data no_var_table;
    set merged(where=(flag_multigene = 0 and flag_no_reads =0 and flag_no_var = 1));
    keep fusion_id;
    run;

proc sort data=no_var_table;
    by fusion_id;
    run;

proc sort data=FRU.results_plus_gov2;
    by fusion_id;
    run;

data merged2;
    merge FRU.results_plus_gov2 (in=in1) no_var_table (in=in2);
    by fusion_id;
    if in2;
    run;

data FRU.dropped_no_var;
   set merged2; 
   keep Fusion_ID CHROM START END symbol_cat Exon_Gene_ID_cat Sequence_loc_cat
   exon_ID_cat Exon_Name_cat FBtrs_per_exon_cat FBgn_cat FBpp_cat FBtr_cat
   Genes_per_fusion Exons_per_fusion FBgns_per_fusion FBpps_per_fusion
   FBtrs_per_fusion Constitutive Common Alternative most_three_prime_exon;
   run;

proc export data=FRU.dropped_no_var outfile='!MCLAB/arbeitman_fru_network/reports_external/fusions_with_no_variance.csv' dbms=CSV replace;
    putnames = yes;
    run;





/* Figure out how many fusions went into the analysis 
    data tmp;
        set FRU.foranalysis_no_multi;
        run;

    proc sort data=tmp nodupkey;
        by fusion_id;
        run; * 45339 fusions;
*/


/* running total 60291 to 45339
    60291 - 2329 = 57962 (minus multigene)
    57962 - 1519 = 56443 (minus no reads)
    56443 - 11104 = 45339 (minus no var)
*/


proc export data=FRU.flag_no_reads outfile='!MCLAB/arbeitman_fru_network/reports_external/flag_no_reads.csv' dbms=CSV replace;
putnames=yes;
run;

