
libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";



data mel_chip_rna_frag_flags_anno3;
merge mel_chip_frag_flags   mel_RNA_frag_detect_flags   mel_RNA_bias    mel_frag_anno;
by featureid;
  /* 216570 obs in all datasets */
RUN;
                                            
proc contents data = mel_chip_rna_frag_flags_anno3 ;
run;

proc freq data = mel_chip_rna_frag_flags_anno3 ;
tables featureType flag_: ;
run;
proc freq data = mel_chip_rna_frag_flags_anno3 ;
tables flag_mel_m_on0_apn * flag_mel_m_on5_apn ;
run;

/*                                            Cumulative
 featureType       Frequency     Percent     Frequency
 ------------------------------------------------------
 3UTR                 21600        9.97         21600
 5UTR                 28479       13.15         50079
 TSS                  22893       10.57         72972
 fragment             87473       40.39        160445
 intergenic           11356        5.24        171801
 intron               44769       20.67        216570
*/

data mel_chip_rna_frag_flags_anno4 ;
set mel_chip_rna_frag_flags_anno3 ;

if featureType="intron" or featureType="intergenic" then flag_coding=0;
	else flag_coding=1;

if featureType="fragment" or featureType="intron" then flag_genic=1;
	else if   featureType="intergenic" then flag_genic=0;

if flag_multigene="" then flag_multigene=-1;    /* make missing = -1 */
    else flag_multigene = flag_multigene ;

if flag_mel_m_on0_apn=1 or flag_mel_f_on0_apn=1 then flag_RNA_detected=1;  
		else flag_RNA_detected=0;
if flag_mel_m_on5_apn=1 or flag_mel_f_on5_apn=1 then flag_RNA_expressed=1;  
		else flag_RNA_expressed=0;

if flag_m_K27_on=1 or flag_f_K27_on=1 then flag_k27_detected=1;
		else flag_k27_detected=0;
if flag_m_K4_on=1 or flag_f_K4_on=1 then flag_k4_detected=1;
		else flag_k4_detected=0;
	
if chrom="X" then xsome="X";
	else if chrom="4" then xsome="4";
	else if chrom="Y" then xsome="Y";
	else if chrom="2L" or chrom="2R" or chrom="3L" or chrom="3R" then xsome="A";
	else xsome="";

if chrom="X" then XCHR=1; 
	else if  xsome="" then xchr="";
	else XCHR=0;
if chrom="4" or chrom="Y" then flag_4Y=1;
	else flag_4Y=0;
if chrom="X" then flag_X=1;
	else flag_X=0;
if chrom="Y" then flag_Y=1;
	else flag_Y=0;

if flag_m_K4_on=1 and flag_f_K4_on=1 then k4_detected="both";
	else if flag_m_K4_on=1 and flag_f_K4_on=0 then k4_detected="male";
	else if flag_m_K4_on=0 and flag_f_K4_on=1 then k4_detected="fem";
	else if flag_m_K4_on=0 and flag_f_K4_on=0 then k4_detected="none";
	else k4_detected="oops";

if flag_m_k27_on=1 and flag_f_k27_on=1 then k27_detected="both";
	else if flag_m_k27_on=1 and flag_f_k27_on=0 then k27_detected="male";
	else if flag_m_k27_on=0 and flag_f_k27_on=1 then k27_detected="fem";
	else if flag_m_k27_on=0 and flag_f_k27_on=0 then k27_detected="none";
	else k27_detected="oops";

if flag_m_K4_on=1 and flag_f_K4_on=1 then k4_sex_ratio="both";
	else if flag_m_K4_on=1 and flag_f_K4_on=0 then k4_sex_ratio="male";
	else if flag_m_K4_on=0 and flag_f_K4_on=1 then k4_sex_ratio="fem";
	else if flag_m_K4_on=0 and flag_f_K4_on=0 then k4_sex_ratio="";
	else k4_sex_ratio="oops";

if flag_m_K27_on=1 and flag_f_K27_on=1 then k27_sex_ratio="both";
	else if flag_m_k27_on=1 and flag_f_k27_on=0 then k27_sex_ratio="male";
	else if flag_m_k27_on=0 and flag_f_k27_on=1 then k27_sex_ratio="fem";
	else if flag_m_k27_on=0 and flag_f_k27_on=0 then k27_sex_ratio="";
	else k27_sex_ratio="oops";

if flag_mel_m_on0_apn=1 and flag_mel_f_on0_apn=0 then rna_detected="male";
	else if flag_mel_m_on0_apn=0 and flag_mel_f_on0_apn=1 then rna_detected="fem";
	else if flag_mel_m_on0_apn=1 and flag_mel_f_on0_apn=1 then rna_detected="both";
	else if flag_mel_m_on0_apn=0 and flag_mel_f_on0_apn=0 then rna_detected="none";
	else  rna_detected="oops";

run ;

proc freq data = mel_chip_rna_frag_flags_anno4 ;
tables 
flag_genic flag_multigene flag_RNA_detected flag_k27_detected flag_k4_detected flag_X flag_Y XCHR k4_detected k27_detected k4_sex_ratio k27_sex_ratio rna_detected
flag_m_bias_RNA_apn_uq_ff flag_f_bias_RNA_apn_uq_ff;
run;


data mel_chip_rna_frag_flags_anno5 ;
set mel_chip_rna_frag_flags_anno4 ;

r_f_m_apn_uq_ff_chk = ratio_avg_f_m_apn_uq_ff ;

if ratio_avg_f_m_apn_uq_ff=1 then r_f_m_apn_uq_ff_chk=.;

if rna_detected="none" then r_f_m_apn_uq_ff_chk=.; /*ensures no bias calls when detection is too low*/

if flag_m_bias_RNA_apn_uq_ff=1 and rna_detected="fem" then r_f_m_apn_uq_ff_chk=.;
if flag_f_bias_RNA_apn_uq_ff=1 and rna_detected="male" then r_f_m_apn_uq_ff_chk=.;

if r_f_m_apn_uq_ff_chk le 0.5 and not missing(r_f_m_apn_uq_ff_chk) then ratio_expressed="male";
    else if r_f_m_apn_uq_ff_chk ge 2.0 then ratio_expressed="fem";
	else  ratio_expressed="unb";
	
if r_f_m_apn_uq_ff_chk=. then ratio_expressed="";

/*new piece of code to identify trend in expression*/
if r_f_m_apn_uq_ff_chk < 1 and not missing(r_f_m_apn_uq_ff_chk) then ratio_trend="male";
    else if r_f_m_apn_uq_ff_chk > 1  then ratio_trend="fem";
	else  ratio_trend="unb";
	
if r_f_m_apn_uq_ff_chk=. then ratio_trend="";

run;

/* check how many frags with no ratio_trend */
proc freq data =  mel_chip_rna_frag_flags_anno5;
tables ratio_trend ;
run;  /*    fem = 64,317
            male = 92,764
            missing = 59,489 */

data check_missing ;
retain ratio_avg_f_m_apn_uq_ff flag_mel_f_allsID_off0_apn flag_mel_m_allsID_off0_apn rna_detected ;
set mel_chip_rna_frag_flags_anno5 ;
where ratio_trend = "";
run;

proc freq data = check_missing ;
where r_f_m_apn_uq_ff_chk = . ;
tables rna_detected / out = checkr ;
tables flag_mel_f_allsID_off0_apn / out = checkf  ;
tables flag_mel_m_allsID_off0_apn / out = checkm  ;
tables flag_mel_f_allsID_off0_apn * flag_mel_m_allsID_off0_apn / out = checkfm  ;
tables ratio_avg_f_m_apn_uq_ff * rna_detected / out = checkdet  ;
run;   /* there are 21026 where r_f_m_apn_uq_ff_chk is missing and */

data check2_missing ;
set check_missing ;
where rna_detected ne "none" ;
run;

proc freq data = check2_missing ;
tables flag_multigene * rna_detected * flag_mel_;
run ;

data check3_missing ;
set check_missing ;
where rna_detected ne "none" and flag_multigene = 0 ;
run;  /* looks ok  - missing are either not expressed or male-female bias is not consistent */



proc sort data=mel_chip_rna_frag_flags_anno5 ;
by fbgn;
run;

proc contents data = mel_chip_rna_frag_flags_anno5 ; run ;

/* do not include anything with ambiguous annotation 
    dropping multigene (includes bidirectional TSS) and fragments with no fb617 annotations */
data mel_chip_rna_frag_flags_anno anno_miss gff_miss oops check_gene;
      /* 216,570 obs                         191,883 obs             12,065 obs */   
merge mel_chip_rna_frag_flags_anno5 (in=in1) fbgn_anno_all (in=in2) fbgn_2go (in=in3);  
by fbgn;
if in1 and in2 then output mel_chip_rna_frag_flags_anno;   /* 198,811 obs */
else if in1 then output  anno_miss;                         /* 17,759 obs */
else if in2 then output gff_miss;                       /* 174,406 pbs */
else if in3 then output check_gene;                     /* 1,183 obs */
else output oops;                                       /* 0 in oops */
run;


PROC EXPORT DATA= WORK.mel_chip_rna_frag_flags_anno 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/data_files/mel_chip_rna_frag_flags_anno.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



data chiprna.mel_chip_rna_frag_flags_anno;
set  mel_chip_rna_frag_flags_anno;
run;

;
