

/*

merge retention time into BFF filtered data

*/

%macro add_rt (type) ;

/* import 'annot' files containing RT */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/text_data/meta_analysis_&type._SLAW_output_0xBFF_PH_annot.tsv"
out = anno_&type.
dbms = tab replace ;
guessingrows = max ;
run ;


/* keep rt and mz columns */
data anno2_&type. ;
retain uniqueID mz rt mean_mz mean_rt max_mz max_rt min_mz min_rt;
set anno_&type. ;
keep uniqueID max_mz max_rt mean_mz mean_rt min_mz min_rt mz rt ;
run ;

/* import bff filtered data */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/bff_filtering/BFF_corrected_&type._PD_poolPD.tsv"
out = raw_&type.
dbms = tab replace ;
guessingrows = max ;
run ;

proc sort data = anno2_&type. ;
by uniqueID ;
proc sort data = raw_&type. ;
by uniqueID ;
run ;

/* merge and keep if in bff filtered data */
data &type._w_rt filtered_&type. ;
merge anno2_&type. (in=in1) raw_&type. (in=in2) ;
by uniqueID ;
if in2 then output &type._w_rt ;
else output filtered_&type. ;
run;

proc export data = &type._w_rt
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/CID/sweet16/SLAW_UGA_Output/bff_filtering/BFF_corrected_&type._with_rt_mz.tsv"
dbms = tab replace ;
run ;
%mend ;

%add_rt (rp_pos) ;
%add_rt (rp_neg) ;
%add_rt (hilic_pos) ;
%add_rt (hilic_neg) ;





