
libname chiprna "!MCLAB/Dros_PB_ChIP/sasdata/chipRnaMs";



ods pdf file= "!MCLAB/Dros_PB_ChIP/CHIP_RNA_ms/results/agreement/mel_chip_agreement_more.pdf" ;


data mel_chip_rna_frag_flags_anno1;
set chiprna.mel_chip_rna_frag_flags_anno;
	where featuretype="fragment" and flag_multigene=0;
	run;


title 'mel fragments no multigene';

proc freq data=chiprna.mel_chip_rna_frag_flags_anno;
tables flag_mel_m_on0_apn*flag_mel_f_on0_apn/agree;
run;

proc freq data=chiprna.mel_chip_rna_frag_flags_anno;
tables featuretype*flag_mel_m_on0_apn*flag_mel_f_on0_apn/agree;
run;


title1 "mel male female by feature type";
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
tables featuretype*(ratio_expressed);
run;

title1 "mel fragments no multigene";
proc freq data=mel_chip_rna_frag_flags_anno1;
	tables xsome*ratio_expressed*rna_detected;
	run;

proc contents data = mel_chip_rna_frag_flags_anno1; run;
/*
proc freq data= mel_chip_rna_frag_flags_anno1;
tables xsome*rna_detected*k4_detected;
tables xsome*ratio_expressed*k4_detected;
tables xsome*sex_limited*k4_detected;
run;

title1 "RNA detected ";
proc freq data= mel_chip_rna_frag_flags_anno1;
where rna_detected ne 'none';
tables xsome*ratio_expressed;
tables xsome*sex_bias;
tables xsome*sex_limited;
run;


title1 "clear sex bias in both chip and rna TSS";
proc freq data= chiprna.mel_chip_rna_frag_flags_anno;
where rna_detected ne 'none' and featureType="TSS";
tables xsome*sex_bias;
tables xsome*sex_limited;
run;
*/
title "";
ods pdf close ;

