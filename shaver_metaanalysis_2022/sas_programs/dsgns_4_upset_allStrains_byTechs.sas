libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

create table containing effect sizes from meta-strain models all strains - omitting N2!!

sig in at least 1 strain (le 0.05)

*/




filename mymacros "!MCLAB/maize_ozone_FINAL/2018/PacBio/sas_programs/macros";
options SASAUTOS=(sasautos mymacros);



data dsgn ;
set sweet16.dsgn_gt_rp_neg_slaw ;
keep strain ;
if strain = "" then delete ;
if strain = "PD1074" then delete ;
run ;

proc sort data = dsgn  nodups  ;
by _all_ ;
run;

/* import mass spec effect sizes */
%macro imp_outer (tech) ;
%macro imp (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_analysis_&tech./MA_FE_rank_byMut_&tech._&strain._SMD_summary.tsv"
out = &strain._&tech. 
dbms = csv replace ;
run;

data a_&tech._&strain.;
set &strain._&tech. ;
rename effect = effect_&strain. ;
rename p_value = pvalue_&strain. ;
if p_value le 0.05 then flag_sig_&strain = 1 ;
else flag_sig_&strain = 0;
keep effect featureID p_value flag_sig_&strain. ;
run ;

proc sort data =a_&tech._&strain.;
by featureID ;
run;

%mend ;

%iterdataset(dataset=dsgn, function=%nrstr(%imp(&strain);)); 
%mend ;

%imp_outer (rp_neg) ;
%imp_outer (rp_pos) ;
%imp_outer (hilic_pos) ;


/* merge all strains by tech */
%macro byTech (tech) ;

data &tech._all ;
merge a_&tech._: ;
by featureID ;
run ;

data sig_&tech. ;
set &tech._all ;
if flag_sig_AUM2073 = 1 or
flag_sig_CB4856 = 1 or
flag_sig_CX11314 = 1 or
flag_sig_DL238 = 1 or
flag_sig_KJ550 = 1 or
flag_sig_N2 = 1 or
flag_sig_RB2011 = 1 or
flag_sig_RB2055 = 1 or
flag_sig_RB2347 = 1 or
flag_sig_RB2550 = 1 or
flag_sig_UGT49 = 1 or
flag_sig_UGT60 = 1 or
flag_sig_VC1265 = 1 or
flag_sig_VC2524 then flag_keep = 1 ;
run;

data upset_&tech ;
set sig_&tech. ;
where flag_keep = 1 ;
keep featureID flag_sig_: ;
run ;


proc export data = upset_&tech
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/upset_plots/MA_byStrain_sig_4_upset_&tech..csv"
dbms = csv replace ;
run;


proc contents data = upset_&tech  out = vars_&tech. noprint;
run;

data upset_dsgn_&tech. ;
set vars_&tech ;
rename name = columnID;
var1 = scan(name, 1, '_') ;
var2 = scan(name, 2, '_') ;
var3 = scan(name, 3, '_') ;
label = compress(var2||'_'||var3) ;
keep name label ;
run ;

proc export data = upset_dsgn_&tech.
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/upset_plots/upset_dsgn_&tech..csv"
dbms = csv replace ;
run;

%mend ;

%byTech (rp_neg);
%byTech (rp_pos);
%byTech (hilic_pos);




/* same for NMR */
/* import effect sizes */
%macro imp_outer (tech) ;
%macro imp (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/nmr/meta_analysis/MA_fixedEffect_NMR_&tech._rank_&strain._&tech._SMD_summary.tsv"
out = &strain._&tech. 
dbms = csv replace ;
run;

data a_&tech._&strain.;
set &strain._&tech. ;
rename effect = effect_&strain. ;
rename p_value = pvalue_&strain. ;
if p_value le 0.05 then flag_sig_&strain = 1 ;
else flag_sig_&strain = 0;
keep effect featureID p_value flag_sig_&strain. ;
run ;

proc sort data =a_&tech._&strain.;
by featureID ;
run;

%mend ;

%iterdataset(dataset=dsgn, function=%nrstr(%imp(&strain);)); 
%mend ;

%imp_outer (cdcl3) ;
%imp_outer (d2o) ;


/* merge all strains by tech */
%macro byTech (tech) ;

data &tech._all ;
merge a_&tech._: ;
by featureID ;
run ;

data sig_&tech. ;
set &tech._all ;
if flag_sig_AUM2073 = 1 or
flag_sig_CB4856 = 1 or
flag_sig_CX11314 = 1 or
flag_sig_DL238 = 1 or
flag_sig_KJ550 = 1 or
flag_sig_RB2011 = 1 or
flag_sig_RB2055 = 1 or
flag_sig_RB2347 = 1 or
flag_sig_RB2550 = 1 or
flag_sig_UGT49 = 1 or
flag_sig_UGT60 = 1 or
flag_sig_VC1265 = 1 or
flag_sig_VC2524 then flag_keep = 1 ;
run;

data upset_&tech ;
set sig_&tech. ;
where flag_keep = 1 ;
keep featureID flag_sig_: ;
run ;

proc export data = upset_&tech
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/nmr/upset_plots/MA_byStrain_sig_4_upset_&tech..csv"
dbms = csv replace ;
run;

proc contents data = upset_&tech  out = vars_&tech. noprint;
run;

data upset_dsgn_&tech. ;
set vars_&tech ;
rename name = columnID;
var1 = scan(name, 1, '_') ;
var2 = scan(name, 2, '_') ;
var3 = scan(name, 3, '_') ;
label = compress(var2||'_'||var3) ;
keep name label ;
run ;

proc export data = upset_dsgn_&tech.
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/nmr/upset_plots/upset_dsgn_&tech..csv"
dbms = csv replace ;
run;

%mend ;

%byTech (cdcl3);
%byTech (d2o);









