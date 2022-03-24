libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*
mass spec pathway NI without N2

rp_neg, rp_pos and hilic_pos

create table containing effect sizes for sig

CB4856, CX11314, DL238
*/

/* import meta-pathway NI no N2 */
%macro impTech (tech) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&tech./MA_FE_rank_byPath_&tech._NInoN2_SMD_summary.tsv"
out = path_&tech.
dbms = csv replace ;
run;

data path2_&tech. ;
set path_&tech. ;
rename p_value = pValue_path ;
keep featureID p_value ;
run;

proc sort data = path2_&tech. ;
by featureID ;
run;

/* import meta-strain - strains in N2 except N2 */
%macro imping (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&tech./MA_FE_rank_byMut_&tech._&strain._SMD_summary.tsv"
out = &strain._&tech.
dbms = csv replace ;
run;

data a_&strain._&tech. ;
set &strain._&tech. ;
rename p_value = pValue_&strain. ;
keep featureID p_value ;
run;

proc sort data = a_&strain._&tech. ;
by featureID ;
run;
%mend ;

%imping (CB4856) ;
%imping (CX11314) ;
%imping (DL238) ;

%mend ;

%impTech (rp_neg) ;
%impTech (rp_pos) ;
%impTech (hilic_pos) ;


%macro forest (tech) ;

data forest_&tech. ;
merge path2_&tech. a_&tech._:;
by featureID ;
run;

data forestPlot_&tech. ;
set forest_&tech. ;
if pvalue_path le 0.05 then flag_path = 1 ; else flag_path = 0 ;
if pvalue_CB4856 le 0.05 then flag_CB4856 = 1 ; else flag_CB4856 =0;
if pvalue_CX11314 le 0.05 then flag_CX11314 = 1 ; else flag_CX11314 =0;
if pvalue_DL238 le 0.05 then flag_DL238 = 1 ; else flag_DL238 = 0;
if flag_path = 1 and flag_CB4856 = 0 and flag_CX11314 = 0 and flag_DL238 = 0 then flag_plot = 1 ;
run ;

data plotit_&tech. ;
set forestPlot_&tech. ;
where flag_plot = 1 ;
run;

data list_&tech. ;
set plotit_&tech. ;
keep featureID ;
run ;

proc export data = list_&tech. 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&tech./list_features_sig_metaPath_NInon2_notSig_metaStrain_&tech..tsv"
dbms = tab replace ;
run;


%mend ;
%forest (rp_neg);
%forest (rp_pos);
%forest (hilic_pos);



