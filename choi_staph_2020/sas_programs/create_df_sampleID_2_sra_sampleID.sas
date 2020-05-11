
libname relapse "/home/ammorse/staph_relapse/sasdata";
libname lmm "/home/ammorse/staph_relapse/sasdata/data_analysis";

/*
create design file to rename sampleID in vcf files to sra_sampleID

Does NOT include failed library and pair (3938 and 3867) 
	78 obs
*/


data relapse.sampleID_2_sra_sampleID;
set relapse.sampleID_st_cc_ref_kitchen_sink ;
if patientnumber = 16 then delete ;
keep sampleID sra_sampleID reference;
run;


proc export data = relapse.sampleID_2_sra_sampleID
outfile = "/home/ammorse/staph_relapse/design_files/sampleID_2_sra_sampleID.csv"
dbms = csv replace ;
run;

proc export data = relapse.sampleID_2_sra_sampleID
outfile = "/home/ammorse/staph_relapse/design_files/sampleID_2_sra_sampleID_noHeader.csv"
dbms = csv replace ;
putnames = no ;
run;
