libname relapse "!MCLAB/staph/seoung_ho_project/sasdata";
libname TB14 "/home/ammorse/TB14/staph_relapse/sasdata";



/*
additional scatter plots
(1)
Yaxis = total# snps relative to ref for A
Xaxis  = total# snps relative to ref for B
(2)
Yaxis = snps in B not in A
Xaxis = snps in A not in B
*/

ods pdf file = "!MCLAB/staph/seoung_ho_project/output/A_B_mismatch_scatter_plot.pdf" ;

title "SNPS in A not B" ;
proc sgplot data = relapse.relapse_snp_cnts_ST_CC_ref_sra ;
where is_a_pair = "YES";
scatter x = genotype_pair_1_0 y = genotype_pair_0_1/ markerattrs=(symbol=CircleFilled size=10px) ;
xaxis label = 'SNPS in A not B (genotype_pair_1_0)' ;
yaxis label = 'SNPS in B not A (genotype_pair_0_1)' ;
run;
ods pdf close ;


data play ;
set relapse.relapse_snp_cnts_ST_CC_ref_sra ;
A_ref = sum(genotype_pair_0_0 + genotype_pair_0_1) ;
B_ref = sum(genotype_pair_0_0 + genotype_pair_1_0) ;
run ;

ods pdf file = "!MCLAB/staph/seoung_ho_project/output/A_ref_vs_B_ref_scatter_plot.pdf" ;

title "SNPs relative to reference (SNPs with no coverage omitted" ;
proc sgplot data = play;
where is_a_pair = "YES";
scatter x = A_ref y = B_ref / markerattrs=(symbol=CircleFilled size=10px) ;
xaxis label = 'SNPS in A relative to ref (genotype_pair_0_0 + genotype_pair_0_1))' ;
yaxis label = 'SNPS in B relative to ref (genotype_pair_0_0 + genotype_pair_1_0))' ;
run;

ods pdf close ;
