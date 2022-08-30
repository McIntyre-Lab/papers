libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

data for upset plots 
    need binary flag for whether feature is sig


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

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_analysis_&type./MA_FE_rank_byMut_&type._&strain._SMD_summary.tsv"
out = &strain._&type.
dbms = csv replace ;
run;

data A_&type._&strain. ;
set &strain._&type. ;
if p_value le 0.05 then flag_&strain._sig = 1 ;
else flag_&strain._sig = 0 ;
if effect ge 0 then flag_&strain._effect_ge_pd1074 = 1;
else flag_&strain._effect_ge_pd1074 = 0;
keep featureID flag_&strain._sig flag_&strain._effect_ge_pd1074;
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

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_analysis_&type./MA_FE_rank_byPath_&type._&path._SMD_summary.tsv"
out = &path._&type. 
dbms = csv replace ;
run;

data A_&type._&path. ;
set &path._&type. ;
if p_value le 0.05 then flag_&path._sig = 1 ;
else flag_&path._sig = 0 ;
if effect ge 0 then flag_&path._effect_ge_pd1074 = 1;
else flag_&path._effect_ge_pd1074 = 0;
keep featureID flag_&path._sig flag_&path._effect_ge_pd1074 ;
run ;

proc sort data = A_&type._&path. ;
by featureID ;
run;

%mend ;

%pvPath (NI);
%pvPath (NInoN2);
%pvPath (TCA);
%pvPath (UGT);
%pvPath (Combo);

%mend ;

%pv_types (rp_neg) ;
%pv_types (rp_pos) ;
%pv_types (hilic_pos) ;




/* merge flags for each tech - pathway combination */

%macro combine (tech) ;


data &tech._TCA_4_upset ;
merge A_&tech._TCA A_&tech._KJ550 A_&tech._RB2347 A_&tech._VC1265 A_&tech._AUM2073 A_&tech._VC2524;
by featureID ;
run;

data &tech._UGT_4_upset ;
merge A_&tech._UGT A_&tech._RB2055 A_&tech._RB2550 A_&tech._RB2011 A_&tech._UGT49 A_&tech._UGT60;
by featureID ;
run;

data &tech._NI_4_upset ;
merge A_&tech._NI A_&tech._CB4856 A_&tech._CX11314 A_&tech._DL238 A_&tech._N2;
by featureID ;
run;

data &tech._NInoN2_4_upset ;
merge A_&tech._NI A_&tech._CB4856 A_&tech._CX11314 A_&tech._DL238 ;
by featureID ;
run;

data &tech._Combo_4_upset ;
merge A_&tech._Combo A_&tech._VC1265 A_&tech._DL238 A_&tech._UGT60;
by featureID ;
run;

%mend ;

%combine (rp_neg);
%combine (rp_pos);
%combine (hilic_pos);



%macro exporting (dataout) ;

proc export data = &dataout. 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/upset_plots/&dataout._plot.csv"
dbms = csv replace ;
run ;
%mend ;

%exporting (rp_neg_TCA_4_upset) ;
%exporting (rp_pos_TCA_4_upset) ;
%exporting (hilic_pos_TCA_4_upset) ;

%exporting (rp_neg_UGT_4_upset) ;
%exporting (rp_pos_UGT_4_upset) ;
%exporting (hilic_pos_UGT_4_upset) ;

%exporting (rp_neg_NI_4_upset) ;
%exporting (rp_pos_NI_4_upset) ;
%exporting (hilic_pos_NI_4_upset) ;

%exporting (rp_neg_NInoN2_4_upset) ;
%exporting (rp_pos_NInoN2_4_upset) ;
%exporting (hilic_pos_NInoN2_4_upset) ;

%exporting (rp_neg_Combo_4_upset) ;
%exporting (rp_pos_Combo_4_upset) ;
%exporting (hilic_pos_Combo_4_upset) ;


%macro dsgns (path) ;

proc contents data = hilic_pos_&path._4_upset out = listing_&path. ; run;

data design_&path. ;
set listing_&path. ;
keep name ;
rename name = columnID ;
if name = "featureID" then delete ;
run ;

data design2_&path. ;
set design_&path. ;
label columnID = "columnID" ;
if find(columnID, "sig") ge 1 then do ;
    var2 = compress(scan(columnID, 2, '_'));
    var3 = compress(scan(columnID, 3, '_'));
    label = compress(var2||'_'||var3);
end ;
if find(columnID, "effect") ge 1 then do ;
    var2 = compress(scan(columnID, 2, '_'));
    var3 = compress(scan(columnID, 3, '_'));
    var4 = compress(scan(columnID, 4, '_'));
    var5 = compress(scan(columnID, 5, '_'));
    label = compress(var2||'_'||var3||'_'||var4||'_'||var5);
end ;
keep columnID label ;
run ;

data dsgn_&path._sig ;
set design2_&path.;
if find(label, "sig") ge 1 ;
run ;

data dsgn_&path._ge_pd1074 ;
set design2_&path.;
if find(label, "effect") ge 1 ;
run ;

proc export data = dsgn_&path._sig
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/upset_plots/upset_plot_design_file_&path._sig.csv"
dbms = csv replace ;
run ;

proc export data = dsgn_&path._ge_pd1074
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/upset_plots/upset_plot_design_file_&path._ge_pd1074.csv"
dbms = csv replace ;
run ;

%mend ;

%dsgns (TCA) ;
%dsgns (UGT) ;
%dsgns (NI) ;
%dsgns (NInoN2) ;
%dsgns (Combo) ;




