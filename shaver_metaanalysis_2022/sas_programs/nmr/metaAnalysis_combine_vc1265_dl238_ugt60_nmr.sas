libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

combine effect size tables for:
    combo pathway (VC1265, DL238, UGT60, PD1074 )
    by strain for VC1265, DL238, and UGT60
    ** 4 tables to merge

NMR DATA
*/





/* %macro import (MA_fixedEffect_pathway_NMR_cdcl3_rank_COMBO_summary, cdcl3, COMBO) ;*/
%macro import (datain, which, name) ;


proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/&datain..tsv"
out = &name._&which.
dbms = csv replace ;
run;

proc sql ;
select cat(name, '=', cats(name, "_&name.")) into :suffix separated by ' ' from 
dictionary.columns where libname = 'WORK' and upcase(memname) = upcase("&name._&which.") ;
quit ;

data a_&which._&name;
set &name._&which. (rename = (&suffix));
run ;

data &which._&name.;
set a_&which._&name; 
rename featureID_&name. = featureID ;
run ;

proc sort data = &which._&name.;
by featureID ;
run;

%mend ;

%import (MA_fixedEffect_pathway_NMR_cdcl3_rank_COMBO_summary, cdcl3, COMBO) ;
%import (MA_fixedEffect_pathway_NMR_d2o_rank_COMBO_summary, d2o, COMBO) ;

%import (MA_fixedEffect_NMR_cdcl3_rank_DL238_cdcl3_SMD_summary, cdcl3, DL238) ;
%import (MA_fixedEffect_NMR_d2o_rank_DL238_d2o_SMD_summary, d2o, DL238) ;

%import (MA_fixedEffect_NMR_cdcl3_rank_UGT60_cdcl3_SMD_summary, cdcl3, UGT60) ;
%import (MA_fixedEffect_NMR_d2o_rank_UGT60_d2o_SMD_summary, d2o, UGT60) ;

%import (MA_fixedEffect_NMR_cdcl3_rank_VC1265_cdcl3_SMD_summary, cdcl3, VC1265) ;
%import (MA_fixedEffect_NMR_d2o_rank_VC1265_d2o_SMD_summary, d2o, VC1265) ;


%macro combine (type) ;

data MA_output_4_combo_&type ;
merge &type._: ;
by featureID ;
run ;

proc export data =  MA_output_4_combo_&type 
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/nmr/meta_analysis/MA_output_4_combo_&type..csv"
dbms = csv replace ;
run ;

%mend ;

%combine(cdcl3) ;
%combine(d2o) ;

