libname sweet16 "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/sasdata" ;


/*

all the annotations to unranked and ranked analysis ready datasets --> needed for metabolomics workbench

format mean_sn max_sn and min_sn so all values in E format 

keep 4 decimal places for mz and 2 for RT

*/




/*  import big anno file */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/text_data/meta_analysis_rp_neg_SLAW_output_0xBFF_PH_annot.tsv"
out = rp_neg_anno 
dbms = tab replace ;
guessingrows = MAX ;
run;

data check ;
retain UniqueID ms2_id min_height mean_sn max_right_on_left_assymetry clique mean_right_on_left_assymetry isotopic_pattern_annot min_mz min_sn max_rt mean_rt_cor 
raw_isotopic_pattern rt min_right_on_left_assymetry min_rt_cor mean_peakwidth mean_mz min_peakwidth max_sn mz max_height max_rt_cor min_rt mean_height 
neutral_mass num_clustered_ms2 min_intensity mean_rt main_peak max_peakwidth max_mz num_detection annotations mean_intensity max_intensity total_detection ;

format max_sn e11.;
format mean_sn e11.;
format min_sn e11.;
set rp_neg_anno ;
run ;

proc contents data = rp_neg_anno order = varnum  ;
run ;

/*  order of variables in original annot -- keep 
UniqueID ms2_id min_height mean_sn max_right_on_left_assymetry clique mean_right_on_left_assymetry isotopic_pattern_annot min_mz min_sn max_rt mean_rt_cor raw_isotopic_pattern rt min_right_on_left_assymetry min_rt_cor mean_peakwidth mean_mz min_peakwidth max_sn mz max_height max_rt_cor min_rt mean_height neutral_mass num_clustered_ms2 min_intensity mean_rt main_peak max_peakwidth max_mz num_detection annotations mean_intensity max_intensity total_detection
*/



%macro wb (datain, tech, techNorm, dataout) ;

/*  import analysis ready dataset */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/analysisReady_datasets/&datain."
out = &techNorm 
dbms = tab replace ;
run;

/*  import big anno file */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/text_data/meta_analysis_&tech._SLAW_output_0xBFF_PH_annot.tsv"
out = &tech._anno 
dbms = tab replace ;
guessingrows = MAX ;
run;

data &tech._anno2 ;
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
set  &tech._anno ;
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

%wb (rp_neg_analysisReady_rank_sbys.tsv, rp_neg, rp_neg_rank, rp_neg_analysisReady_rank_sbys_wAnno_full.tsv ) ;
%wb (rp_neg_analysisReady_sbys.tsv, rp_neg, rp_neg, rp_neg_analysisReady_sbys_wAnno_full.tsv) ;

%wb (rp_pos_analysisReady_rank_sbys.tsv, rp_pos, rp_pos_rank, rp_pos_analysisReady_rank_sbys_wAnno_full.tsv ) ;
%wb (rp_pos_analysisReady_sbys.tsv, rp_pos, rp_pos, rp_pos_analysisReady_sbys_wAnno_full.tsv) ;

%wb (hilic_pos_analysisReady_rank_sbys.tsv, hilic_pos, hilic_pos_rank, hilic_pos_analysisReady_rank_sbys_wAnno_full.tsv ) ;
%wb (hilic_pos_analysisReady_sbys.tsv, hilic_pos, hilic_pos, hilic_pos_analysisReady_sbys_wAnno_full.tsv) ;


