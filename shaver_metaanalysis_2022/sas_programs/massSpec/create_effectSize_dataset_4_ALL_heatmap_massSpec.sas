

/*
mass spec effect sizes

datasets for heatmap of ALL pathway with all mutants

*/


%macro pv_types (type) ;

%macro pvalues (strain) ;
/* import MA output files */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_analysis_&type./MA_FE_rank_byMut_&type._&strain._SMD_summary.tsv"
out = &strain._&type.
dbms = csv replace ;
run;

/* flag if strain is sig */
data A_&type._&strain. ;
set &strain._&type. ;
rename effect = effect_&strain. ;
keep featureID effect ;
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
rename effect = effect_&path._path ;
keep featureID effect ;
run ;

proc sort data = A_&type._&path. ;
by featureID ;
run;

%mend ;
%pvPath (ALL);

%mend ;

%pv_types (rp_neg) ;
%pv_types (rp_pos) ;
%pv_types (hilic_pos) ;



%macro combine (tech) ;

data MA_effectSizes_&tech._ALL ;
merge A_&tech._ALL A_&tech._AUM2073 A_&tech._CB4856 A_&tech._CX11314 A_&tech._DL238  A_&tech._KJ550 A_&tech._N2  A_&tech._RB2011 A_&tech._RB2055 A_&tech._RB2347 A_&tech._RB2550 A_&tech._UGT49 A_&tech._UGT60 A_&tech._VC1265 A_&tech._VC2524 ;;
by featureID ;
run;

proc export data = MA_effectSizes_&tech._ALL 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/shaver_metaanalysis_2022/massSpec/meta_plotting/MA_effectSizes_&tech._ALL.csv"
dbms = csv replace ;
run ;

%mend ;

%combine (rp_neg);
%combine (rp_pos);
%combine (hilic_pos);





