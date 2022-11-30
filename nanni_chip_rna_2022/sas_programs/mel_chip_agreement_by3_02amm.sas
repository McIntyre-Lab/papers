
libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";




proc freq data=chiprna.mel_chip_rna_frag_flags_anno;
tables featuretype;
run;

ods pdf file= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_chip_agreement.pdf" ;

title "MEL male female by feature type";
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
tables featuretype*(flag_k4_detected flag_k27_detected);
run;

/*green figure*/
title "MEL k27 irrelevant k4 detected";
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
where featuretype ne "intergenic" and k4_detected ne "none";
tables featuretype*k4_detected;
run;
title "MEL only k4 detected";
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
where featuretype ne "intergenic" and k4_detected ne "none" and flag_k27_detected=0;
tables featuretype*k4_detected;
run;

proc sort data= chiprna.mel_chip_rna_frag_flags_anno;
by featuretype;
run;

proc contents data = chiprna.mel_chip_rna_frag_flags_anno ;
run ;

proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
by featuretype;
tables flag_m_K4_on * flag_f_K4_on/agree;
run;

proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
by featuretype;
where flag_k27_detected=0;
tables flag_m_K4_on*flag_f_K4_on/agree;
run;

/* k4_on does not exist and k27_on does not exist 
proc contents data = chiprna.mel_chip_rna_frag_flags_anno; run;
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
*where k27_on="none" and k4_on ne "none";
tables featuretype*k4_on/chisq;
run;
*/

proc freq data=chiprna.mel_chip_rna_frag_flags_anno;
tables k4_sex_ratio;
tables featuretype*k4_sex_ratio;
run;

/* k4_on does not exist and k27_on does not exist 
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
where k4_on ne "none" and (featureType="3UTR" or featureType="fragment");
tables featuretype*k4_on*k27_on/cmh;
run;

proc sort data=chiprna.mel_chip_rna_frag_flags_anno;
by xsome;
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
by xsome;
tables featuretype*k4_on*k27_on;
run;
*/

proc sort data=chiprna.mel_chip_rna_frag_flags_anno;
by featuretype xsome;
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
by featuretype xsome;
tables flag_m_k4_on*flag_f_k4_on/agree out=mel_ft_k4;
 output out=mel_agree_ft_k4 kappa mcnem;
run;


proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
by featuretype xsome;
tables flag_m_k27_on*flag_f_k27_on/agree out=mel_ft_k27 ;
 output out=mel_agree_ft_k27  kappa mcnem;
run;



title "mel f_v_m k4 no k27";
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
by featuretype xsome;
where k27_detected="none";
tables flag_m_k4_on*flag_f_k4_on/agree out=mel_ft_k4_nok27;
output out=mel_agree_ft_k4_nok27 kappa mcnem;
run;


proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
by featuretype xsome;
where k4_detected="none";
tables flag_m_k27_on*flag_f_k27_on/agree out=mel_ft_k27_nok4;
output out=mel_agree_ft_k27_nok4 kappa mcnem;
run;


title 'mel';

proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
where flag_genic=1 and (xsome="X" or xsome="A");
tables xsome*k4_detected*k27_detected/cmh;
run;

proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
where flag_genic=1 and (xsome="X" or xsome="A");
tables xsome*k4_detected/chisq;
tables xsome*k27_detected/chisq;
tables xsome*k4_sex_ratio/chisq;
tables xsome*k27_sex_ratio/chisq;
run;

title '';

ods pdf close ;


PROC EXPORT DATA= WORK.mel_agree_ft_k4 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_agree_ft_k4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.mel_agree_ft_k27
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_agree_ft_k27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.mel_ft_k4 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_ft_k4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.mel_ft_k27
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_ft_k27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


PROC EXPORT DATA= WORK.mel_agree_ft_k4_nok27  
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_agree_ft_k4_nok27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.mel_agree_ft_k27_nok4
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_agree_ft_k27_nok4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.mel_ft_k4_nok27  
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_ft_k4_nok27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.mel_ft_k27_nok4
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_ft_k27_nok4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


