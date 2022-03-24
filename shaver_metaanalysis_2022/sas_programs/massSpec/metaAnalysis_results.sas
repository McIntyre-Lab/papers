libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

Effect size table by pathway

    if group == 'TCA':
        muts = ["KJ550", "RB2347", "VC1265", "AUM2073", "VC2524"]
    elif group == 'UGT':
        muts = ["RB2055", "RB2550", "RB2011", 'UGT49', 'UGT60']
    elif group == 'NI':
        muts = ["CB4856", "CX11314", "DL238", "N2"]


*/


data list ;
set sweet16.dsgn_gt_rp_neg_slaw ;
keep strain ;
run ;

proc sort data = list  nodups ; ;
by _all_ ;
run;

proc print data = list ; run;



%macro pv_types (type) ;

%macro pvalues (strain, path) ;
/* import MA output files */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&type./MA_FE_rank_byMut_&type._&strain._SMD_summary.tsv"
out = &strain._&type.
dbms = csv replace ;
run;

/* flag if strain is sig */
data A_&type._&strain. ;
set &strain._&type. ;
if p_value le 0.05 then flag_sig05_&strain. = 1 ;
else flag_sig05_&strain. = 0 ;

if effect > 0 then direction_&strain = "strain";
else if effect le 0 then direction_&strain = "pd1074";

rename effect = effect_&strain. ;
keep featureID flag_sig05_&strain. effect direction_&strain ;
run ;

proc sort data = A_&type._&strain. ;
by featureID ;
run;

%mend ;

%pvalues (AUM2073);
%pvalues (CB4856);
%pvalues (CX11314);
%pvalues (DL238);
%pvalues (KJ550);
%pvalues (N2);
%pvalues (RB2011);
%pvalues (RB2055);
%pvalues (RB2347);
%pvalues (RB2550);
%pvalues (UGT49);
%pvalues (UGT60);
%pvalues (VC1265);
%pvalues (VC2524);


%mend ;

%pv_types (rp_neg) ;
%pv_types (rp_pos) ;
%pv_types (hilic_pos) ;




%macro pv_types (type) ;

%macro pvPath (path) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&type./MA_FE_rank_byPath_&type._&path._SMD_summary.tsv"
out = &path._&type. 
dbms = csv replace ;
run;

data A_&type._&path. ;
set &path._&type. ;
if p_value le 0.05 then flag_sig05_&path._path = 1 ;
else flag_sig05_&path._path  = 0 ;
rename effect = effect_&path._path ;
keep featureID flag_sig05_&path._path effect ;
run ;

proc sort data = A_&type._&path. ;
by featureID ;
run;

%mend ;

%pvPath (NI);
%pvPath (TCA);
%pvPath (UGT);

%mend ;

%pv_types (rp_neg) ;
%pv_types (rp_pos) ;
%pv_types (hilic_pos) ;





/* merge flags for each tech - pathway combination */

%macro combine (tech) ;

data MA_output_&tech._TCA ;
merge A_&tech._TCA A_&tech._KJ550 A_&tech._RB2347 A_&tech._VC1265 A_&tech._AUM2073 A_&tech._VC2524;
by featureID ;
run;

data MA_output_&tech._UGT ;
merge A_&tech._UGT A_&tech._RB2055 A_&tech._RB2550 A_&tech._RB2011 A_&tech._UGT49 A_&tech._UGT60;
by featureID ;
run;

data MA_output_&tech._NI ;
merge A_&tech._NI A_&tech._CB4856 A_&tech._CX11314 A_&tech._DL238 A_&tech._N2;
by featureID ;
run;
%mend ;

%combine (rp_neg);
%combine (rp_pos);
%combine (hilic_pos);


/* combination flags */

%macro flags (tech) ;

data MA_output_&tech._NI_flags ;
set  MA_output_&tech._NI ;
length direction_pd_vs_strain $24. ;
sum_flag_sig05_strain = sum(flag_sig05_CB4856 + flag_sig05_CX11314 + flag_sig05_DL238 + flag_sig05_N2) ; 

if direction_CB4856 = "strain" and direction_CX11314 = "strain" and direction_DL238 = "strain" and direction_N2 = "strain" then direction_pd_vs_strain = "all strain direction" ;
else if direction_CB4856 = "pd1074" and direction_CX11314 = "pd1074" and direction_DL238 = "pd1074" and direction_N2 = "pd1074" then direction_pd_vs_strain = "all pd1074 direction" ;
else direction_pd_vs_strain = "mix strain & pd1074" ;
run ;


data MA_output_&tech._UGT_flags ;
set  MA_output_&tech._UGT ;
length direction_pd_vs_strain $24. ;
sum_flag_sig05_strain = sum(flag_sig05_RB2055 + flag_sig05_RB2550 + flag_sig05_RB2011 + flag_sig05_UGT49 + flag_sig05_UGT60) ; 

if direction_RB2055 = "strain" and direction_RB2550 = "strain" and direction_RB2011 = "strain" and direction_UGT49 = "strain" and direction_UGT60 = "strain" then direction_pd_vs_strain = "all strain direction" ;
else if direction_RB2055 = "pd1074" and direction_RB2550 = "pd1074" and direction_RB2011 = "pd1074" and direction_UGT49 = "pd1074" and direction_UGT60 = "pd1074" then direction_pd_vs_strain = "all pd1074 direction" ;
else direction_pd_vs_strain = "mix strain & pd1074" ;
run ;


data MA_output_&tech._TCA_flags ;
set  MA_output_&tech._TCA ;
length direction_pd_vs_strain $24. ;
sum_flag_sig05_strain = sum(flag_sig05_KJ550 + flag_sig05_RB2347 + flag_sig05_VC1265 + flag_sig05_AUM2073 + flag_sig05_VC2524) ; 

if direction_KJ550 = "strain" and direction_RB2347 = "strain" and direction_VC1265 = "strain" and direction_AUM2073 = "strain" and direction_VC2524 = "strain" then direction_pd_vs_strain = "all strain direction" ;
else if direction_KJ550 = "pd1074" and direction_RB2347 = "pd1074" and direction_VC1265 = "pd1074" and direction_AUM2073 = "pd1074" and direction_VC2524= "pd1074" then direction_pd_vs_strain = "all pd1074 direction" ;
else direction_pd_vs_strain = "mix" ;
run ;

%mend ;

%flags (rp_neg);
%flags (rp_pos);
%flags (hilic_pos);


%macro exporting (dataout) ;

proc export data = &dataout. 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/comparing_models/&dataout..csv"
dbms = csv replace ;
run ;
%mend ;

%exporting (MA_output_rp_neg_NI_flags) ;
%exporting (MA_output_rp_neg_TCA_flags) ;
%exporting (MA_output_rp_neg_UGT_flags) ;

%exporting (MA_output_rp_pos_NI_flags) ;
%exporting (MA_output_rp_pos_TCA_flags) ;
%exporting (MA_output_rp_pos_UGT_flags) ;

%exporting (MA_output_hilic_pos_NI_flags) ;
%exporting (MA_output_hilic_pos_TCA_flags) ;
%exporting (MA_output_hilic_pos_UGT_flags) ;

