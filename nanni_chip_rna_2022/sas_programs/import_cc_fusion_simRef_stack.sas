libname ortho "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/sas_data/RNAseq_ortho";
libname df  "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/svn_lmm_dros_head_data/ethanol_srna/sas_data";


filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



/*
create stacked cc dataset 
*/



/* list of samples to loop over format - mel_301_m_noEtoh_rep2 */
data dsgn2 ;
set df.design_ethanol_srna_rnaseq ;
keep species genotype sex arbitrary_rep treatment;
run ;

proc sort data = dsgn2 nodups ;
by _all_ ;
run;

data dsgn ;
retain datain sample ;
set dsgn2 ;
sample = compress(species||'_'||genotype||'_'||sex||'_'||treatment||'_rep'||arbitrary_rep) ;
if species = "sim" then end = "combined" ;
else if species = "mel" then end = "mel2sim" ;
datain = compress('cvrg_cnts_'||sample||'_'||end) ;
keep datain sample ;
run;  /* 48 samples */


proc sort data = dsgn nodups ;
by _all_ ;
run; /* uniq */

data test_list ;
set dsgn ;
where sample = "sim_11_f_etoh_rep1" ;
run ;


%macro importing (datain, sample) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq/Coverage_counts/fusion_cvrg_cnts_simRef/&datain..csv"
out = in_&sample.
dbms = csv replace ;
guessingrows = MAX;
run;   /* 64,619 features */

data in2_&sample ;
retain sampleID ;
set in_&sample ;
sampleID = "&sample" ;
run;

%mend ;

%iterdataset(dataset=dsgn, function=%nrstr(%importing(&datain, &sample);)); 


data cvrg_fusion_stack;
length sampleID $21. ;
set in2_: ;
run ;  /* 64378 obs per in2 file */

data ortho.cvrg_fusion_simRef_stack ;
set cvrg_fusion_stack ;
run;

proc freq data = ortho.cvrg_fusion_simRef_stack ;
tables sampleID  /out = cnt_samples;
run ;

