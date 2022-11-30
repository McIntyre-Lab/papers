
libname ortho "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/sas_data/RNAseq_ortho";
libname drosRNA "!MCLAB/Dros_CHIP_RNA_ms/sas_data/RNAseq";
libname df  "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/svn_lmm_dros_head_data/ethanol_srna/sas_data";


filename mymacros "!MCLAB/maize_ozone_final/2014/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/*
on - offs for fusion coverage counts of mel to sim and sim to sim
    flag a fusion as 'on' if in both species for each sex (flag_sim_F_on = 1)
    line and trt are extra reps

    flag on if 50% of reps are on at apn > 5 --> flag_&species._&geno._&sex._&trt._on5 = 1

into norm (log uq) if flag_sim_f_on = 1 and flag_mel_f_on = 1, etc

*/

data dsgn2 ;
set df.design_ethanol_srna_rnaseq ;
keep species genotype sex arbitrary_rep treatment;
run ;

proc sort data = dsgn2 nodups ;
by _all_ ;
run;

data design ;
retain sampleID ;
length sampleID $21. ;
set dsgn2 ;
sampleID = compress(species||'_'||genotype||'_'||sex||'_'||treatment||'_rep'||arbitrary_rep) ;
sample = compress(species||'_'||sex) ;
rep = compress('rep'||arbitrary_rep);
drop arbitrary_rep ;
run;  /* 48 samples */

proc sort data = ortho.cvrg_fusion_simRef_stack ;
by sampleID ;
proc sort data  = design;
by sampleID ;
run;

data cnts_df oops  ;
merge ortho.cvrg_fusion_simRef_stack  (in=in1) design (in=in2) ;
by sampleID ;
if in1 and in2 then output cnts_df ;
else output oops ;
run;  /* 0 obs in oops */

data cnt_flag ;
set cnts_df;
if apn > 0 then apn_on0 = 1; else apn_on0 = 0;
if apn > 5 then apn_on5 = 1; else apn_on5 = 0;
run;


%macro onCalls (apn) ;

proc sort data=cnt_flag;
  by featureID species sex  ;
  run;

proc means data = cnt_flag noprint ;
    by featureID species sex ;
    var apn_on&apn ;
    output out = sample_on&apn mean=sample_percent_on ;
    run ;

proc sort data = sample_on&apn ;
    by featureID ;
    run;

data sample2_on&apn ;
retain featureID sample ;
set  sample_on&apn ;
/* fusion is on if expressed in in 50% of the reps */
if sample_percent_on > 0.5 then flag_fusion_on&apn = 1 ;
    else flag_fusion_on&apn = 0 ;
sample = compress(species||'_'||sex);
run ;


proc transpose data = sample2_on&apn out = sample_sbys_on&apn prefix = flag_ suffix= _on&apn;
    by featureID ;
    id sample ;
    var flag_fusion_on&apn ;
    run;

data ortho.onCall_sample_APN&apn._simRef ;
set sample_sbys_on&apn;
if flag_mel_F_on&apn. = 1 and flag_sim_F_on&apn. = 1 then flag_F_on&apn = 1; else flag_F_on&apn = 0 ;
if flag_mel_M_on&apn. = 1 and flag_sim_M_on&apn. = 1 then flag_M_on&apn = 1; else flag_M_on&apn = 0 ;
drop _name_ ;
run;

%mend  ;
%onCalls (0) ;
%onCalls (5) ;


proc freq data = ortho.onCall_sample_APN0_simRef ;
tables flag_F_on0 * flag_M_on0 ;
run;  /* 64,378 fusions */
/*  
 flag_F_on0     flag_M_on0

 Frequency|
 Percent  |
 Row Pct  |
 Col Pct  |       0|       1|  Total
 ---------+--------+--------+
        0 |  13285 |   1291 |  14576
          |  20.64 |   2.01 |  22.64
          |  91.14 |   8.86 |
          |  94.82 |   2.56 |
 ---------+--------+--------+
        1 |    726 |  49076 |  49802
          |   1.13 |  76.23 |  77.36
          |   1.46 |  98.54 |
          |   5.18 |  97.44 |
 ---------+--------+--------+
 Total       14011    50367    64378
             21.76    78.24   100.00


*/
