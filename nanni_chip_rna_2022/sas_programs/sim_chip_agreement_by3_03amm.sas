


libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";

ods pdf file= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_chip_agreement.pdf" ;

title "sim male female by feature type";
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
tables featuretype*(flag_k4_detected flag_k27_detected);
run;
title "SIM k27 irrelevant  k4 detected";
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
where featuretype ne "intergenic" and k4_detected ne "none";
tables featuretype*k4_detected;
run;

title "SIM only k4 detected";
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
where featuretype ne "intergenic" and k4_detected ne "none" and flag_k27_detected=0;
tables featuretype*k4_detected;
run;

proc sort data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype;
run;

proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype;
tables flag_m_K4_on*flag_f_K4_on/agree;
run;


proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype;
where flag_k27_detected=0;
tables flag_m_K4_on*flag_f_K4_on/agree;
run;

/* k4_on does not exist and k27_on does not exist 
proc sort data=chiprna.sim_chip_rna_frag_flags_anno;
by xsome;
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by xsome;
tables featuretype*k4_on*k27_on;
run;
*/

proc sort data=chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
tables flag_m_k4_on*flag_f_k4_on/agree;
 output out=sim_agree_k4 kappa mcnem;
run;


proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
tables flag_m_k27_on*flag_f_k27_on/agree;
 output out=sim_agree_k27 kappa mcnem;
run;



proc sort data=chiprna.sim_chip_rna_frag_flags_anno;
by  xsome;
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by xsome;
tables flag_m_k4_on*flag_f_k4_on/agree;
* output out=agree_ft_k4 kappa mcnem;
run;


proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by xsome;
tables flag_m_k27_on*flag_f_k27_on/agree;
*output out=agree_ft_k27 kappa mcnem;
run;


proc sort data=chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
tables flag_m_k4_on*flag_f_k4_on/agree out=sim_ft_k4;
 output out=sim_agree_ft_k4 kappa mcnem;
run;


proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
tables flag_m_k27_on*flag_f_k27_on/agree out=sim_ft_k27;
 output out=sim_agree_ft_k27 kappa mcnem;
run;


proc sort data=chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;

title "sim f_v_m k4 no k27";
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
where k27_detected="none";
tables flag_m_k4_on*flag_f_k4_on/agree out=sim_ft_k4_nok27;
output out=sim_agree_ft_k4_nok27 kappa mcnem;
run;


proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
by featuretype xsome;
where k4_detected="none";
tables flag_m_k27_on*flag_f_k27_on/agree out=sim_ft_k27_nok4;
output out=sim_agree_ft_k27_nok4 kappa mcnem;
run;


title 'sim';
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
where flag_genic=1 and (xsome="X" or xsome="A");
tables xsome*k4_detected*k27_detected/cmh;
run;


title 'sim';
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
where flag_genic=1 and (xsome="X" or xsome="A");
tables xsome*k4_detected/chisq;
tables xsome*k27_detected/chisq;
tables xsome*k4_sex_ratio/chisq;
tables xsome*k27_sex_ratio/chisq;
run;
title "";
ods pdf close ;



PROC EXPORT DATA= WORK.sim_agree_k4 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_agree_k4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


PROC EXPORT DATA= WORK.sim_agree_k27 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_agree_k27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.sim_agree_ft_k4 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_agree_ft_k4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.sim_agree_ft_k27 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_agree_ft_k27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.sim_ft_k4 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_ft_k4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.sim_ft_k27 
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_ft_k27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;





PROC EXPORT DATA= WORK.sim_agree_ft_k4_nok27  
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_agree_ft_k4_nok27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.sim_agree_ft_k27_nok4
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_agree_ft_k27_nok4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


PROC EXPORT DATA= WORK.sim_ft_k4_nok27  
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_ft_k4_nok27.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



PROC EXPORT DATA= WORK.sim_ft_k27_nok4
            OUTFILE= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_ft_k27_nok4.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

