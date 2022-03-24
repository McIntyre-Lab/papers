libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

add mz and RT to unranked and ranked analysis ready datasets --> needed for metabolomics workbench

keep 4 decimal places for mz and 2 for RT

*/


%macro wb (datain, tech, techNorm, dataout) ;

/*  import analysis ready dataset */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/analysisReady_datasets/&datain."
out = &techNorm 
dbms = tab replace ;
run;

/*  import anno file with mz and rt */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/filtered_datasets_processed/BFF_corrected_&tech._with_rt_mz_solv_filt_annot.tsv"
out = &tech._anno 
dbms = tab replace ;
run;

proc contents data = rp_neg_anno ; run;

data &tech._anno2 ;
retain uniqueID mz RT ;
format mz 12.4 ;
format rt 12.2 ;
set  &tech._anno ;
label mz = "m/z" ;
label rt = "RT" ;
keep uniqueID mz rt ;
run ;

proc sort data = &techNorm ;
by uniqueID ;
proc sort data = &tech._anno2 ;
by uniqueID ;
run ;

data &techNorm._wb ;
merge &tech._anno2 (in=in1) &techNorm. (in=in2) ;
by uniqueID ;
if in2 ;
run ;
proc export data = &techNorm._wb
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/analysisReady_datasets/&dataout."
dbms = tab replace ;
run ;

%mend ;

%wb (rp_neg_analysisReady_rank_sbys.tsv, rp_neg, rp_neg_rank, rp_neg_analysisReady_rank_sbys_wAnno.tsv ) ;
%wb (rp_neg_analysisReady_sbys.tsv, rp_neg, rp_neg, rp_neg_analysisReady_sbys_wAnno.tsv) ;

%wb (rp_pos_analysisReady_rank_sbys.tsv, rp_pos, rp_pos_rank, rp_pos_analysisReady_rank_sbys_wAnno.tsv ) ;
%wb (rp_pos_analysisReady_sbys.tsv, rp_pos, rp_pos, rp_pos_analysisReady_sbys_wAnno.tsv) ;

%wb (hilic_pos_analysisReady_rank_sbys.tsv, hilic_pos, hilic_pos_rank, hilic_pos_analysisReady_rank_sbys_wAnno.tsv ) ;
%wb (hilic_pos_analysisReady_sbys.tsv, hilic_pos, hilic_pos, hilic_pos_analysisReady_sbys_wAnno.tsv) ;


