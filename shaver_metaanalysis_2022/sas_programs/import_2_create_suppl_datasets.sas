libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;

filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);

/*

create supplemental files 
    (1) import analysis ready datasets, import MA results, create sig flags, import all annotations, sirius output, etc

*/

data list ;
set sweet16.dsgn_gt_rp_neg_slaw ;
keep strain ;
if strain = "" then delete ;
if strain = "PD1074" then delete ;
run ;

proc sort data = list  nodups  ;
by _all_ ;
run;


/* import mass spec effect sizes and pvalues*/
%macro imp_outer (tech) ;
%macro imp (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&tech./MA_FE_rank_byMut_&tech._&strain._SMD_summary.tsv"
out = &strain._&tech. 
dbms = csv replace ;
run;

data ma_&tech._&strain.;
set &strain._&tech. ;
rename effect = effect_&strain. ;
rename p_value = pval_&strain. ;
if p_value le 0.05 then flag_sig_&strain = 1 ;
else flag_sig_&strain = 0;
keep featureID effect p_value flag_sig_&strain ;
run ;

proc sort data =ma_&tech._&strain.;
by featureID ;
run;

%mend ;

%iterdataset(dataset=list, function=%nrstr(%imp(&strain);)); 
%mend ;

%imp_outer (rp_neg) ;
%imp_outer (rp_pos) ;
%imp_outer (hilic_pos) ;


/* run for pathway */
%macro imp_outer2 (tech) ;
%macro imp2 (path) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&tech./MA_FE_rank_byPath_&tech._&path._SMD_summary.tsv"
out = &path._&tech. 
dbms = csv replace ;
run;

data ma_&tech._&path.;
set &path._&tech. ;
rename effect = effect_&path. ;
rename p_value = pval_&path. ;
if p_value le 0.05 then flag_sig_&path = 1 ;
else flag_sig_&path = 0;
keep featureID effect p_value flag_sig_&path ;
run ;

proc sort data = ma_&tech._&path.;
by featureID ;
run;

%mend ;
%imp2 (UGT) ;
%imp2 (TCA) ;
%imp2 (NI) ;
%mend ;

%imp_outer2 (rp_neg) ;
%imp_outer2 (rp_pos) ;
%imp_outer2 (hilic_pos) ;

/* import NMR effect sizes and pvalues */
%macro imp_outer (tech) ;
%macro imp (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/MA_fixedEffect_NMR_&tech._rank_&strain._&tech._SMD_summary.tsv"
out = &strain._&tech. 
dbms = csv replace ;
run;

data ma_&tech._&strain.;
set &strain._&tech. ;
rename effect = effect_&strain. ;
rename p_value = pval_&strain. ;
if p_value le 0.05 then flag_sig_&strain = 1 ;
else flag_sig_&strain = 0;
keep featureID effect p_value flag_sig_&strain ;

proc sort data =ma_&tech._&strain.;
by featureID ;
run;

%mend ;

%iterdataset(dataset=list, function=%nrstr(%imp(&strain);)); 
%mend ;

%imp_outer (d2o) ;
%imp_outer (cdcl3) ;

/* import NMR for pathway */
%macro imp_outer2 (tech) ;
%macro imp2 (path) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/MA_fixedEffect_pathway_NMR_&tech._rank_&path._summary.tsv"
out = &path._&tech. 
dbms = csv replace ;
run;

data ma_&tech._&path.;
set &path._&tech. ;
rename effect = effect_&path. ;
rename p_value = pval_&path. ;
if p_value le 0.05 then flag_sig_&path = 1 ;
else flag_sig_&path = 0;
keep featureID effect p_value flag_sig_&path ;

proc sort data =ma_&tech._&path.;
by featureID ;
run;

%mend ;
%imp2 (UGT) ;
%imp2 (TCA) ;
%imp2 (NI) ;
%mend ;

%imp_outer2 (d2o) ;
%imp_outer2 (cdcl3) ;



/* import PH and annotations for mass spec */


%macro wb (datain, tech, techNorm) ;

/*  import analysis ready dataset */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/analysisReady_datasets/&datain."
out = ph2_&techNorm 
dbms = tab replace ;
run;

data ph_&techNorm ;
set ph2_&techNorm ;
rename uniqueID = featureID ;
run;

/*  import big anno file */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/text_data/meta_analysis_&tech._SLAW_output_0xBFF_PH_annot.tsv"
out = an2_&tech._anno 
dbms = tab replace ;
guessingrows = MAX ;
run;

data an_&tech._anno ;
retain UniqueID ms2_id min_height mean_sn max_right_on_left_assymetry clique mean_right_on_left_assymetry isotopic_pattern_annot min_mz min_sn max_rt mean_rt_cor 
raw_isotopic_pattern rt min_right_on_left_assymetry min_rt_cor mean_peakwidth mean_mz min_peakwidth max_sn mz max_height max_rt_cor min_rt mean_height 
neutral_mass num_clustered_ms2 min_intensity mean_rt main_peak max_peakwidth max_mz num_detection annotations mean_intensity max_intensity total_detection ;

format max_sn e11.;
format mean_sn e11.;
format min_sn e11.;
format mz 12.4 ;
format max_mz  12.4 ;
format min_mz 12.4 ;
format mean_mz  12.4 ;
format rt 12.2 ;
format max_rt 12.2 ;
format min_rt  12.2 ;
format mean_rt  12.2 ;
rename uniqueID = featureID ;
set  an2_&tech._anno ;
run ;

proc sort data = ph_&techNorm ;
by featureID ;
proc sort data = an_&tech._anno ;
by featureID ;
run ;

%mend ;


%wb (rp_neg_analysisReady_sbys.tsv, rp_neg, rp_neg) ;
%wb (rp_pos_analysisReady_sbys.tsv, rp_pos, rp_pos) ;
%wb (hilic_pos_analysisReady_sbys.tsv, hilic_pos, hilic_pos) ;




/* import nmr sbys data */
%macro imp_ph2 (tech ) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/ready_&tech._PH_sbys.tsv"
out = ph2_&tech 
dbms = tab replace ;
guessingrows = MAX ;
run;

data ph_&tech. ;
set ph2_&tech. ;
rename ppm = featureID ;
run ;

proc sort data = ph_&tech. ;
by featureID ;
run;

%mend ;

%imp_ph2 (cdcl3) ;
%imp_ph2 (d2o) ;


/* import SIRIUS data (4 mass spec) */
%macro imp_sirius (tech );

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/MA_SIRIUS_output/SIRIUS_output_mergedIDs/SIRIUS_MA_&tech._mergedID.csv"
out = sr2_&tech 
dbms = csv replace ;
guessingrows = MAX ;
run;

data sr_&tech. ;
set sr2_&tech.;
rename uniqueID = featureID ;
drop ms2_id num_strains_sig ;  /* drop ms2_id - already in annotations */
run;

proc sort data = sr_&tech. ;
by featureID ;
run ;

%mend ;

%imp_sirius (rp_pos) ; /* 1050 obs */
%imp_sirius (rp_neg) ; /* 211 obs */
%imp_sirius (hilic_pos) ; /* 166 obs */



