
libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";



ods pdf file= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/sim_chip_agreement_more.pdf" ;



data sim_chip_rna_frag_flags_anno1;
set chiprna.sim_chip_rna_frag_flags_anno;
	where featuretype="fragment" and flag_multigene=0;
	run;


title 'sim fragments no multigene';

proc freq data=chiprna.sim_chip_rna_frag_flags_anno;
tables flag_sim_m_on0_apn*flag_sim_f_on0_apn/agree;
run;

proc freq data=chiprna.sim_chip_rna_frag_flags_anno;
tables featuretype*flag_sim_m_on0_apn*flag_sim_f_on0_apn/agree;
run;


title1 "sim male female by feature type";
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
tables featuretype*(ratio_expressed);
run;

title1 "sim fragments no multigene";
proc freq data=sim_chip_rna_frag_flags_anno1;
	tables xsome*ratio_expressed*rna_detected;
	run;

/*
proc freq data= sim_chip_rna_frag_flags_anno1;
tables xsome*rna_detected*k4_detected;
tables xsome*ratio_expressed*k4_detected;
tables xsome*sex_limited*k4_detected;
run;

title1 "sim RNA detected ";
proc freq data= sim_chip_rna_frag_flags_anno1;
where rna_detected ne 'none';
tables xsome*ratio_expressed;
tables xsome*sex_bias;
tables xsome*sex_limited;
run;


title1 "sim clear sex bias in both chip and rna TSS";
proc freq data= chiprna.sim_chip_rna_frag_flags_anno;
where rna_detected ne 'none' and featureType="TSS";
tables xsome*sex_bias;
tables xsome*sex_limited;
run;
*/
title "";

ods pdf close ;
