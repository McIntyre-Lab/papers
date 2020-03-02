
libname tappas "!MCLAB/maize_ozone_FINAL/2018/PacBio/sasdata/tappas";

title "";
/*
looking at data 

create pie chart of regualted genes across 5 genotypes

prep data for hierarchical clustering (heatmap)
*/



/* for genes regulated in all 5 genotypes, what is pattern? */
data tappas_genes_all ;
retain flag_in_all ;
set tappas.tappas_results_genes ;
if DEA_result_B73 ne "" and DEA_result_C123 ne "" and DEA_result_Hp301 ne "" and DEA_result_Mo17 ne "" and DEA_result_NC338 ne "" then flag_in_all = 1 ;
else flag_in_all =0 ;
if regulated_B73 ne '' then DEA_result_B73 = regulated_B73 ; else DEA_result_B73 = DEA_result_B73 ;
if regulated_C123 ne '' then DEA_result_C123 = regulated_C123 ; else DEA_result_C123 = DEA_result_C123 ;
if regulated_Hp301 ne '' then DEA_result_Hp301 = regulated_Hp301 ; else DEA_result_Hp301 = DEA_result_Hp301 ;
if regulated_Mo17 ne '' then DEA_result_Mo17 = regulated_Mo17 ; else DEA_result_Mo17 = DEA_result_Mo17 ;
if regulated_NC338 ne '' then DEA_result_NC338 = regulated_NC338 ; else DEA_result_NC338 = DEA_result_NC338 ;
run ;

proc freq data = tappas_genes_all;
tables DEA_result_: ;
run ;
proc freq data = tappas_genes_all ;
where flag_in_all = 1;
tables DEA_result_B73 * DEA_result_C123 * DEA_result_Hp301 * DEA_result_Mo17 * DEA_result_NC338  / out = dea_gene_freqs ;
run ;



/* for regulated genes, what is the overlap between the 5 genotypes? */

data gene_cnts2 ;
set tappas.tappas_results_genes ;
where regulated_B73  ne "" or regulated_Mo17  ne "" or regulated_C123  ne "" or regulated_Hp301  ne "" or regulated_NC338  ne "" ; 
keep gene regulated_: ;
run;

data gene_cnts ;
set gene_cnts2 ;
if regulated_B73 ne ""  then flag_on_B73 = 1 ;     else flag_on_B73 = 0 ;
if regulated_Mo17 ne "" then flag_on_Mo17 = 1 ;   else flag_on_Mo17 = 0 ;
if regulated_C123 ne "" then flag_on_C123 = 1 ;   else flag_on_C123 = 0 ;
if regulated_Hp301 ne "" then flag_on_Hp301 = 1 ; else flag_on_Hp301 = 0 ;
if regulated_NC338 ne "" then flag_on_NC338 = 1 ; else flag_on_NC338 = 0 ;
counts = (flag_on_B73 + flag_on_Mo17 + flag_on_C123 + flag_on_Hp301 + flag_on_NC338) ;
drop regulated_:;
run ;


proc freq data = gene_cnts ;
tables flag_on_B73 * flag_on_C123 * flag_on_Mo17 * flag_on_Hp301 * flag_on_NC338 / out = gene_venn ;
run;
    
data gene_venn_cnts ;
set gene_venn ;
*drop percent ;
sum_on = (flag_on_B73 + flag_on_C123 + flag_on_Mo17 + flag_on_Hp301 + flag_on_NC338) ;
run;

proc sort data = gene_venn_cnts ;
by sum_on ;
run;

proc means data = gene_venn_cnts sum ;
var count ;
by sum_on ;
output out = geno_sums sum = sum ;
run;
    /*  regulated in 1 genotype = 1934
                     2 genotype = 1129
                     3 genotype = 861
                     4 genotype = 917
                     5 genotype = 81
                    total      4922 genes at 5% FDR     
*/  
                                


data geno_sum_venn ;
set geno_sums ;
if sum_on = 1 then genoNum = '1 Genotype';
else if sum_on = 2 then genoNum = '2 Genotype';
else if sum_on = 3 then genoNum = '3 Genotype';
else if sum_on = 4 then genoNum = '4 Genotype';
else genoNum = '5 Genotype';
run;


goptions reset = all;

/* legend */
legend1 label = none
    position= (left middle)
    offset= (4,)
    across=1
    value=(color=black)
    shape=bar(4, 1.5) ;

title "Overlap in Regulated Genes (5% FDR, TappAS) Between the 5 Genotypes" ;
proc gchart data = geno_sum_venn ;
pie genoNum / sumvar = sum
    other = 0 
    ascending
    legend = legend1
    value = outside 
    coutline = black
    plabel=(font='Albany AMT/bold' height=1) 
    noheading ;
run ;
quit;



/* add geno counts to 1 genotype */
goptions reset = all;

data test ;
set gene_cnts ;
if counts = 1 then genoNum = '1 genotype' ; 
    else if counts =2 then  genoNum = '2 genotype' ; 
    else if counts =3 then  genoNum = '3 genotype' ; 
    else if counts =4 then  genoNum = '4 genotype' ; 
    else  genoNum = '5 genotype' ; 
run;

proc transpose data = test out = tall ;
by 

/* legend */
legend1 label = none
    position= (left middle)
    offset= (4,)
    across=1
    value=(color=black)
    shape=bar(4, 1.5) ;


ods pdf file = "!MCLAB/maize_ozone_FINAL/pacbio_paper/figs/overlap_test.pdf" ;

title "Overlap in Regulated Genes (5%, TappAS) Between the 5 Genotypes" ;
proc gchart data = gene_cnts ;
donut genoNum / sumvar = sum
    subgroup = genotype 
    other = 0 
    ascending
    legend = legend1
    value = outside 
    coutline = black
    plabel=(font='Albany AMT/bold' height=1) 
    noheading ;
run ;
quit;



