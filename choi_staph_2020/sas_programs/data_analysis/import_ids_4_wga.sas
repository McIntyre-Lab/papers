
libname wgs "Z:\SHARE\McIntyre_Lab\staph\seoung_ho_project\sasdata";

data wgs;
set wgs.sampleid_st_cc_ref_4_sra;
rename sampleid=id
pairnum=sequence_pair;
run;
