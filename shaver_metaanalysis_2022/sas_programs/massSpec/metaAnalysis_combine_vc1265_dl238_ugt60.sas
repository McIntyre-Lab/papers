libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

combine effect size tables for:
    combo pathway (VC1265, DL238, UGT60, PD1074 )
    by strain for VC1265, DL238, and UGT60
    ** 4 tables to merge

*/


data list ;
set sweet16.dsgn_gt_rp_neg_slaw ;
keep strain ;
run ;

proc sort data = list  nodups ; ;
by _all_ ;
run;

proc print data = list ; run;



%macro imp_outer (type) ;
%macro imp (which, which2) ;

proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/meta_analysis_&type./MA_FE_rank_&which._&type._&which2._SMD_summary.tsv"
out = &which2._&type. 
dbms = csv replace ;
run;

proc sql ;
select cat(name, '=', cats(name, "_&which2.")) into :suffix separated by ' ' from 
dictionary.columns where libname = 'WORK' and upcase(memname) = upcase("&which2._&type.") ;
quit ;

data a_&which2._&type.;
set &which2._&type. (rename = (&suffix));
run ;

data &type._&which2.;
set a_&which2._&type.; 
rename featureID_&which2. = featureID ;
run ;

proc sort data = &type._&which2.;
by featureID ;
run;

%mend ;
%imp (byPath, COMBO) ;
%imp (byMut, VC1265) ;
%imp (byMut, DL238) ;
%imp (byMut, UGT60) ;

%mend ;

%imp_outer (rp_neg) ;
%imp_outer (rp_pos) ;
%imp_outer (hilic_pos) ;
%mend;


%macro combine (type) ;

data MA_output_4_combo_&type ;
merge &type._: ;
by featureID ;
run ;

proc export data =  MA_output_4_combo_&type 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/comparing_models/MA_output_4_combo_&type..csv"
dbms = csv replace ;
run ;

%mend ;

%combine(rp_neg) ;
%combine(rp_pos) ;
%combine(hilic_pos) ;



