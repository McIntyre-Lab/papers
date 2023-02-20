
libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";



data sim_chip_rna_frag_flags_anno3;
merge sim_chip_frag_flags sim_RNA_frag_detect_flags sim_RNA_bias sim_frag_anno;
by featureid;

if num_fbgn=1 then flag_single_xscrpt=1;
	else if num_fbgn=. then flag_single_xscrpt=.;
	else if num_fbgn>1 then flag_single_xscrpt=0;
	else flag_single_xscrpt=-1;
run;

proc freq data = sim_chip_rna_frag_flags_anno3 ;
tables featureType chrom ;
run;
/*                                           Cumulative
featureType       Frequency     Percent     Frequency
-----------------------------------------------------
3UTR                 16231        7.91         16231
5UTR                 25081       12.22         41312
TSS                  21069       10.27         62381
fragment             79405       38.70        141786
intergenic           16174        7.88        157960
intron               47236       23.02        205196
*/

data sim_chip_rna_frag_flags_anno4 ;
set sim_chip_rna_frag_flags_anno3 ;

if featureType="intron" or featureType="intergenic" then flag_coding=0;
	else flag_coding=1;

if featureType="fragment" or featureType="intron"  then flag_genic=1;
	else if  featureType="intergenic" then flag_genic=0;
     

if flag_multigene="" then flag_multigene=-1;
    else flag_multigene = flag_multigene ;

if flag_sim_m_on0_apn=1 or flag_sim_f_on0_apn=1 then flag_RNA_detected=1;
		else flag_RNA_detected=0;
if flag_sim_m_on5_apn=1 or flag_sim_f_on5_apn=1 then flag_RNA_expressed=1;
		else flag_RNA_expressed=0;

if flag_m_K27_on=1 or flag_f_K27_on=1 then flag_k27_detected=1;
		else flag_k27_detected=0;
if flag_m_K4_on=1 or flag_f_K4_on=1 then flag_k4_detected=1;
		else flag_k4_detected=0;
	
if chrom="Scf_X" then xsome="X";
	else if chrom="Scf_4" then xsome="4";
	else if chrom="Scf_Y" then xsome="Y";
	else if chrom="Scf_2L" or chrom="Scf_2R" or chrom="Scf_3L" or chrom="Scf_3R" then xsome="A";
	else xsome="";

if chrom="Scf_X" then XCHR=1; 
	else if  xsome="" then xchr="";
	else XCHR=0;
if chrom="Scf_4" or chrom="Y" then flag_4Y=1;
	else flag_4Y=0;
if chrom="Scf_X" then flag_X=1;
	else flag_X=0;
if chrom="Scf_Y" then flag_Y=1;
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

if flag_sim_m_on0_apn=1 and flag_sim_f_on0_apn=0 then rna_detected="male";
	else if flag_sim_m_on0_apn=0 and flag_sim_f_on0_apn=1 then rna_detected="fem";
	else if flag_sim_m_on0_apn=1 and flag_sim_f_on0_apn=1 then rna_detected="both";
	else if flag_sim_m_on0_apn=0 and flag_sim_f_on0_apn=0 then rna_detected="none";
	else  rna_detected="oops";


r_f_m_apn_uq_ff_chk=ratio_avg_f_m_apn_uq_ff;

if ratio_avg_f_m_apn_uq_ff=1 then r_f_m_apn_uq_ff_chk=.;
if rna_detected="none" then r_f_m_apn_uq_ff_chk=.; /*ensures no bias calls when detection is too low*/

if flag_m_bias_RNA_apn_uq_ff=1 and rna_detected="fem" then r_f_m_apn_uq_ff_chk=.;
if flag_f_bias_RNA_apn_uq_ff=1 and rna_detected="male" then r_f_m_apn_uq_ff_chk=.;

if r_f_m_apn_uq_ff_chk le 0.5 and not missing(r_f_m_apn_uq_ff_chk)  then ratio_expressed="male";
	else if r_f_m_apn_uq_ff_chk ge 2.0 then ratio_expressed="fem";
	else  ratio_expressed="unb";

if r_f_m_apn_uq_ff_chk=. then ratio_expressed="";

run;

proc freq data = sim_chip_rna_frag_flags_anno4 ;
 tables ratio_expressed ratio_expressed15 sex_bias flag_X flag_Y  ;
run;

    /*  nothing on Y, 2274 on 4Y
    */
        

/* check how many frags with no ratio_trend */
proc freq data =  sim_chip_rna_frag_flags_anno4;
tables ratio_trend ;
run;  /*    fem = 61572
            male = 82825
            missing = 60799 */
proc sort data=sim_chip_rna_frag_flags_anno4 ;
by fbgn;
run;

/* do not include anything with ambiguous annotation 
    dropping multigene (includes bidirectional TSS) and fragments with no fb617 annotations */
data sim_chip_rna_frag_flags_anno   anno_miss   gff_miss    oops;
/*     205,196 obs                             191,883 obs   */   
merge sim_chip_rna_frag_flags_anno4 (in=in1) fbgn_anno_all (in=in2);
by fbgn;
if in1 and in2 then output sim_chip_rna_frag_flags_anno;    /*  186,186 obs */
else if in1 then output anno_miss ; /*    19,010 obs */
else if in2 then output gff_miss;   /*  176,617 obs */
else output oops;           /* 0 in oops */
run;



PROC EXPORT DATA= WORK.sim_chip_rna_frag_flags_anno 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/data_files/sim_chip_rna_frag_flags_anno.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


data chiprna.sim_chip_rna_frag_flags_anno;
set  sim_chip_rna_frag_flags_anno;
run;



