libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*
NMR pathway NI without N2
cdcl3 and d2o

create table containing effect sizes 

CB4856, CX11314, DL238
*/

/* import meta-pathway NI no N2 */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/MA_fixedEffect_pathway_NMR_d2o_rank_NInoN2_summary.tsv"
out = path
dbms = csv replace ;
run;

data path2 ;
set path ;
rename p_value = pValue_path ;
keep featureID p_value ;
run;

proc sort data = path2 ;
by featureID ;
run;

/* import meta-strain - strains in N2 except N2 */
%macro imping (strain) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/MA_fixedEffect_NMR_d2o_rank_&strain._d2o_SMD_summary.tsv"
out = &strain.
dbms = csv replace ;
run;

data a_&strain ;
set &strain ;
rename p_value = pValue_&strain. ;
keep featureID p_value ;
run;

proc sort data = a_&strain. ;
by featureID ;
run;
%mend ;

%imping (CB4856) ;
%imping (CX11314) ;
%imping (DL238) ;

data forest ;
merge path2 a_:;
by featureID ;
run;

data forestPlot ;
set forest ;
if pvalue_path le 0.05 then flag_path = 1 ; else flag_path = 0 ;
if pvalue_CB4856 le 0.05 then flag_CB4856 = 1 ; else flag_CB4856 =0;
if pvalue_CX11314 le 0.05 then flag_CX11314 = 1 ; else flag_CX11314 =0;
if pvalue_DL238 le 0.05 then flag_DL238 = 1 ; else flag_DL238 = 0;
if flag_path = 1 and flag_CB4856 = 0 and flag_CX11314 = 0 and flag_DL238 = 0 then flag_plot = 1 ;
run ;

data plotit ;
set forestPlot ;
where flag_plot = 1 ;
run;

data list ;
set plotit ;
keep featureID ;
run ;

proc export data = list 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/list_features_sig_metaPath_NInon2_not_sig_metaStrain.tsv"
dbms = tab replace ;
run;



