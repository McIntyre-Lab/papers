


/* import gene summary */
PROC IMPORT OUT= mel_gene
            DATAFILE= "!MCLAB/sex_specific_splicing/zenodo/summary_files/gene_summary_dmel.csv"
      DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

data mel_gene2;
set mel_gene;
rename geneid=dmel_geneid;
run;

/* import master list */
PROC IMPORT OUT= list_of_lists
            DATAFILE= "!MCLAB/useful_dmel_data/gene_lists/master_list_of_lists/master_gene_list_14amm.csv"
     DBMS=CSV REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
     GUESSINGROWS=max; 
RUN;

proc contents data = list_of_lists ; run;


data list2;
set list_of_lists;
dmel_geneid=primary_fbgn;
run;

proc sort data=list2;
by dmel_geneid;
run;

/* merge together - keep if in gene summary */
data mel_list_compare;
merge mel_gene2 (in=in1)   list2;
by dmel_geneid;
if in1;
run;

proc freq data=mel_list_compare;
tables nanni23_dmel_gene_: ;
run;


data mel_list_compare1;
set mel_list_compare;

if telonis09_opp_exon_pattern  ="yes" then telonis09_opp_ep=1;
if telonis09_sig_sexbyprobe ="yes" then telonis09_sexbyp=1;


if nanni23_dmel_gene_ratio2_ttest= "" then nanni23_ratio=.;
	else if nanni23_dmel_gene_ratio2_ttest= "unbiased" then nanni23_ratio=0;
	else  nanni23_ratio=1;
if nanni23_dmel_gene_trend_ttest= "" then nanni23_trend=.;
	else if nanni23_dmel_gene_trend_ttest= "unbiased" then nanni23_trend=0;
	else  nanni23_trend=1;
run;
	
proc contents data = mel_list_compare1 ;
run ;

data mel_list_compare2 ;
set mel_list_compare1 ;	
if nanni23_dmel_gene_ratio2_ttest= "unbiased" then nanni23_ratio_ttest="unbiased";
if nanni23_dmel_gene_ratio2_ttest= "male_and_female" then nanni23_ratio_ttest="B";
if nanni23_dmel_gene_ratio2_ttest= "male" then nanni23_ratio_ttest="M";
if nanni23_dmel_gene_ratio2_ttest= "female" then nanni23_ratio_ttest="F";

if nanni23_dmel_gene_trend_ttest= "unbiased" then nanni23_trend_ttest="unbiased";
if nanni23_dmel_gene_trend_ttest= "male_and_female" then nanni23_trend_ttest="B";
if nanni23_dmel_gene_trend_ttest= "male" then nanni23_trend_ttest="M";
if nanni23_dmel_gene_trend_ttest= "female" then nanni23_trend_ttest="F";

run;
	

proc contents data = mel_list_compare2 ; run ;

%macro compare_list (flag);

proc freq data=mel_list_compare2;
where sexbias ne 'not_evaluated';
tables sexbias*(&flag)/out=count1;
run;

data count2;
set count1;
where &flag=1;
drop percent &flag;
run;

proc transpose data=count2 out=flip_count2;
id sexbias;
run;

data &flag;
set flip_count2;
list="&flag";
drop _name_ _label_;

proc datasets library=work;
	delete flip_count2 count2 count1;
run;

%mend compare_list;

%compare_list(alekseyenko2020_MSL_binding);
%compare_list(arbeitman16_AH_BER_sexBiased);
%compare_list(arbeitman16_AH_CS_sexBiased);
%compare_list(arbeitman16_F_BER_dsx_nullSig);
%compare_list(arbeitman16_F_CS_dsx_nullSig);
%compare_list(arbeitman16_M_BER_dsx_nullSig);
%compare_list(arbeitman16_M_CS_dsx_nullSig);
%compare_list(arbeitman16_dsx_regulated);
%compare_list(chang11_ds_tra);
%compare_list(chang11_f_biased);
%compare_list(chang11_m_biased);
%compare_list(chang11_sex_biased);
%compare_list(chang11_tra_bs);
%compare_list(dalton13_fruMA_overexp_F_induced);
%compare_list(dalton13_fruMA_overexp_F_repress);
%compare_list(dalton13_fruMA_overexp_M_induced);
%compare_list(dalton13_fruMA_overexp_M_repress);
%compare_list(dalton13_fruMB_overexp_F_induced);
%compare_list(dalton13_fruMB_overexp_F_repress);
%compare_list(dalton13_fruMB_overexp_M_induced);
%compare_list(dalton13_fruMB_overexp_M_repress);
%compare_list(dalton13_fruMC_overexp_F_induced);
%compare_list(dalton13_fruMC_overexp_F_repress);
%compare_list(dalton13_fruMC_overexp_M_induced);
%compare_list(dalton13_fruMC_overexp_M_repress);
%compare_list(dalton13_fruM_F_regulation);
%compare_list(dalton13_fruM_M_regulation);
%compare_list(dalton13_fruM_overexp_F_induced);
%compare_list(dalton13_fruM_overexp_F_repress);
%compare_list(dalton13_fruM_overexp_M_induced);
%compare_list(dalton13_fruM_overexp_M_repress);
%compare_list(fear15_genes_added_to_sexdet);
%compare_list(fear16_sig_matedOnly);
%compare_list(fear16_sig_virginOnly);
%compare_list(flyAtl2_Adult_F_Head_enrich_gt);
%compare_list(flyAtl2_Adult_F_Head_enrich_gt2);
%compare_list(flyAtl2_Adult_F_Head_enrich_gt5);
%compare_list(flyAtl2_Adult_F_Head_xpress);
%compare_list(flyAtl2_Adult_M_Head_enrich_gt);
%compare_list(flyAtl2_Adult_M_Head_enrich_gt2);
%compare_list(flyAtl2_Adult_M_Head_enrich_gt5);
%compare_list(flyAtl2_Adult_M_Head_xpress);
%compare_list(flydivas_12spp_pos12);
%compare_list(flydivas_12spp_pos78);
%compare_list(flydivas_12spp_pos88a);
%compare_list(flydivas_melgroup_pos12);
%compare_list(flydivas_melgroup_pos78);
%compare_list(flydivas_melgroup_pos88a);
%compare_list(flydivas_melsubgroup_pos12);
%compare_list(flydivas_melsubgroup_pos78);
%compare_list(flydivas_melsubgroup_pos88a);
%compare_list(garud_2015_h12_top50);
%compare_list(goldman07_ds_tra);
%compare_list(goldman07_ds_tra_not_ds_dsx_fru);
%compare_list(goldman07_sex_bias);
%compare_list(graze12_AI_all_mel_bias);
%compare_list(graze12_AI_all_sim_bias);
%compare_list(graze12_AI_any_mel_bias);
%compare_list(graze12_AI_any_sim_bias);
%compare_list(graze12_AI_partial_mel_bias);
%compare_list(graze12_AI_partial_sim_bias);
%compare_list(graze14_sim_f_bias);
%compare_list(graze14_sim_m_bias);
%compare_list(hartmann11_sexBias_all);
%compare_list(hartmann11_sexBias_body);
%compare_list(hartmann11_sexBias_head);
%compare_list(hatje13);
%compare_list(hudry16_sde_v_adult_midgut);
%compare_list(hudry19_M_v_adult_midgut);
%compare_list(innocenti_morrow_antagonistic);
%compare_list(innocenti_morrow_f_fitness);
%compare_list(innocenti_morrow_m_fitness);
%compare_list(innocenti_morrow_mf_fitness);
%compare_list(kim2003_b52_splicing);
%compare_list(kopp08_olfactory);
%compare_list(luo11_dsx_bs);
%compare_list(mcintyre06_sex_bias);

%compare_list(nanni23_ratio);
%compare_list(nanni23_trend);

%compare_list(newell16_f_bias_TRAP);
%compare_list(newell16_f_bias_input);
%compare_list(newell16_m_bias_TRAP);
%compare_list(newell16_m_bias_input);
%compare_list(newell16_mixed_bias_TRAP);
%compare_list(newell16_mixed_bias_input);
%compare_list(palm21_1012dayFFruP1CnH3K27me3);
%compare_list(palm21_1012dayFFruP1CnH3K4me3);
%compare_list(palm21_1012dayFFruP1PoolH3K27me3);
%compare_list(palm21_1012dayFFruP1PoolH3K4me3);
%compare_list(palm21_1012dayMFruP1CnH3K27me3);
%compare_list(palm21_1012dayMFruP1CnH3K4me3);
%compare_list(palm21_1012dayMFruP1PoolH3K27me3);
%compare_list(palm21_1012dayMFruP1PoolH3K4me3);
%compare_list(palm21_1_dayFElavCnH3K27me3);
%compare_list(palm21_1_dayFElavCnH3K4me3);
%compare_list(palm21_1_dayFElavPoolH3K27me3);
%compare_list(palm21_1_dayFElavPoolH3K4me3);
%compare_list(palm21_1_dayFFruP1CnH3K27me3);
%compare_list(palm21_1_dayFFruP1CnH3K4me3);
%compare_list(palm21_1_dayFFruP1PoolH3K27me3);
%compare_list(palm21_1_dayFFruP1PoolH3K4me3);
%compare_list(palm21_1_dayMElavCnH3K27me3);
%compare_list(palm21_1_dayMElavCnH3K4me3);
%compare_list(palm21_1_dayMElavPoolH3K27me3);
%compare_list(palm21_1_dayMElavPoolH3K4me3);
%compare_list(palm21_1_dayMFruP1CnH3K27me3);
%compare_list(palm21_1_dayMFruP1CnH3K4me3);
%compare_list(palm21_1_dayMFruP1PoolH3K27me3);
%compare_list(palm21_1_dayMFruP1PoolH3K4me3);
%compare_list(redfly_TF);
%compare_list(redfly_TF_target);
%compare_list(ruzika19_num_ant_missVars);
%compare_list(ruzika19_num_ant_nonMissVars);
%compare_list(ruzika19_num_nonAnt_missVars);
%compare_list(ruzika19_num_nonAnt_nonMissVars);
%compare_list(sex_det_pathway);
%compare_list(telonis09_opp_ep);
/* %compare_list(telonis09_opp_exon_pattern); */
%compare_list(telonis09_sexbyp);
/* %compare_list(telonis09_sig_sexbyprobe); */
/*%compare_list(wayne07_annotation_ID);
%compare_list(wayne07_o_pattern_female);
%compare_list(wayne07_o_pattern_male);
%compare_list(wayne07_o_probf_gca_female);
%compare_list(wayne07_o_probf_gca_male); */
%compare_list(wayne07_o_sig_female);
%compare_list(wayne07_o_sig_male);
/* %compare_list(wayne07_on_group);*/ 
%compare_list(zhao17_adaptEvol);
%compare_list(zhao17_adaptEvol_posSel);

data list_compare;
length list $32;
set 
alekseyenko2020_MSL_binding
arbeitman16_AH_BER_sexBiased
arbeitman16_AH_CS_sexBiased
arbeitman16_F_BER_dsx_nullSig
arbeitman16_F_CS_dsx_nullSig
arbeitman16_M_BER_dsx_nullSig
arbeitman16_M_CS_dsx_nullSig
arbeitman16_dsx_regulated
chang11_ds_tra
chang11_f_biased
chang11_m_biased
chang11_sex_biased
chang11_tra_bs
dalton13_fruMA_overexp_F_induced
dalton13_fruMA_overexp_F_repress
dalton13_fruMA_overexp_M_induced
dalton13_fruMA_overexp_M_repress
dalton13_fruMB_overexp_F_induced
dalton13_fruMB_overexp_F_repress
dalton13_fruMB_overexp_M_induced
dalton13_fruMB_overexp_M_repress
dalton13_fruMC_overexp_F_induced
dalton13_fruMC_overexp_F_repress
dalton13_fruMC_overexp_M_induced
dalton13_fruMC_overexp_M_repress
dalton13_fruM_F_regulation
dalton13_fruM_M_regulation
dalton13_fruM_overexp_F_induced
dalton13_fruM_overexp_F_repress
dalton13_fruM_overexp_M_induced
dalton13_fruM_overexp_M_repress
fear15_genes_added_to_sexdet
fear16_sig_matedOnly
fear16_sig_virginOnly
flyAtl2_Adult_F_Head_enrich_gt
flyAtl2_Adult_F_Head_enrich_gt2
flyAtl2_Adult_F_Head_enrich_gt5
flyAtl2_Adult_F_Head_xpress
flyAtl2_Adult_M_Head_enrich_gt
flyAtl2_Adult_M_Head_enrich_gt2
flyAtl2_Adult_M_Head_enrich_gt5
flyAtl2_Adult_M_Head_xpress
flydivas_12spp_pos12
flydivas_12spp_pos78
flydivas_12spp_pos88a
flydivas_melgroup_pos12
flydivas_melgroup_pos78
flydivas_melgroup_pos88a
flydivas_melsubgroup_pos12
flydivas_melsubgroup_pos78
flydivas_melsubgroup_pos88a
garud_2015_h12_top50
goldman07_ds_tra
goldman07_ds_tra_not_ds_dsx_fru
goldman07_sex_bias
graze12_AI_all_mel_bias
graze12_AI_all_sim_bias
graze12_AI_any_mel_bias
graze12_AI_any_sim_bias
graze12_AI_partial_mel_bias
graze12_AI_partial_sim_bias
graze14_sim_f_bias
graze14_sim_m_bias
hartmann11_sexBias_all
hartmann11_sexBias_body
hartmann11_sexBias_head
hatje13
hudry16_sde_v_adult_midgut
hudry19_M_v_adult_midgut
innocenti_morrow_antagonistic
innocenti_morrow_f_fitness
innocenti_morrow_m_fitness
innocenti_morrow_mf_fitness
kim2003_b52_splicing
kopp08_olfactory
luo11_dsx_bs
mcintyre06_sex_bias

nanni23_ratio
nanni23_trend

newell16_f_bias_TRAP
newell16_f_bias_input
newell16_m_bias_TRAP
newell16_m_bias_input
newell16_mixed_bias_TRAP
newell16_mixed_bias_input
palm21_1012dayFFruP1CnH3K27me3
palm21_1012dayFFruP1CnH3K4me3
palm21_1012dayFFruP1PoolH3K27me3
palm21_1012dayFFruP1PoolH3K4me3
palm21_1012dayMFruP1CnH3K27me3
palm21_1012dayMFruP1CnH3K4me3
palm21_1012dayMFruP1PoolH3K27me3
palm21_1012dayMFruP1PoolH3K4me3
palm21_1_dayFElavCnH3K27me3
palm21_1_dayFElavCnH3K4me3
palm21_1_dayFElavPoolH3K27me3
palm21_1_dayFElavPoolH3K4me3
palm21_1_dayFFruP1CnH3K27me3
palm21_1_dayFFruP1CnH3K4me3
palm21_1_dayFFruP1PoolH3K27me3
palm21_1_dayFFruP1PoolH3K4me3
palm21_1_dayMElavCnH3K27me3
palm21_1_dayMElavCnH3K4me3
palm21_1_dayMElavPoolH3K27me3
palm21_1_dayMElavPoolH3K4me3
palm21_1_dayMFruP1CnH3K27me3
palm21_1_dayMFruP1CnH3K4me3
palm21_1_dayMFruP1PoolH3K27me3
palm21_1_dayMFruP1PoolH3K4me3
redfly_TF
redfly_TF_target
ruzika19_num_ant_missVars
ruzika19_num_ant_nonMissVars
ruzika19_num_nonAnt_missVars
ruzika19_num_nonAnt_nonMissVars
sex_det_pathway
telonis09_opp_ep
/* telonis09_opp_exon_pattern */
telonis09_sexbyp
/* telonis09_sig_sexbyprobe */
/*wayne07_annotation_ID
wayne07_o_pattern_female
wayne07_o_pattern_male
wayne07_o_probf_gca_female
wayne07_o_probf_gca_male */
wayne07_o_sig_female
wayne07_o_sig_male
/* wayne07_on_group*/ 
zhao17_adaptEvol
zhao17_adaptEvol_posSel
;
run;

proc sort data = list_compare ;
by list ;
run ;

data list_compare_percent;
set list_compare;
total=B+M+F+unbiased;
per_b=B/total;
per_f=F/total;
per_m=m/total;
per_unbiased=unbiased/total;

expected_b=.037*total;
expected_F=.1765*total;
expected_M=.1765*total;
expected_unb=.61*total;


chisq_overall=(((b-expected_b)*(b-expected_b))/expected_b)+ (((f-expected_f)*(f-expected_f))/expected_f)+(((m-expected_M)*(m-expected_M))/expected_m)+(((unbiased-expected_unb)*(unbiased-expected_unb))/expected_unb);
pvalue_chis_overall=1-probchi(chisq_overall,3);
chisq_both=(((b-expected_b)*(b-expected_b))/expected_b)+(((unbiased-expected_unb)*(unbiased-expected_unb))/expected_unb);
pvalue_chis_both=1-probchi(chisq_both,1);
run;


proc freq data=mel_list_compare1;
where sexbias ne 'not_evaluated';
tables sexbias;
run;

data find ;
set list_compare_percent ;

where list ? "hao" ;
run;

PROC EXPORT DATA= list_compare_percent
            OUTFILE= "/nfshome/ammorse/mnt/ufgi.ahc.ufl.edu-ufgi$/SHARE/McIntyre_Lab/sex_specific_splicing/Tables/list_compare_percent_03amm.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;












