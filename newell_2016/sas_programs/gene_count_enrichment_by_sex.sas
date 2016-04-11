
*count unique male bias genes;
proc freq data=all_results noprint;
where male_enrich=1 and female_enrich=0;
tables symbol_cat/out=count_male_enrich;
run;

*count unqiue female biased genes;
proc freq data=all_results noprint;
where female_enrich=1 and male_enrich=0;
tables symbol_cat/out=count_female_enrich;
run;
proc freq data=all_results noprint;
where female_enrich=1 and male_enrich=1;
tables symbol_cat/out=count_enrich_both;
run;

data gene_enrich_by_sex;
merge count_male_enrich (in=in1)  count_female_enrich(in=in2)
count_enrich_both (in=in3);
by symbol_cat;
if in1 and in2 and in3 then gene_compare_enrichment_bysex="male_female_both";
else if in1 and in2 then  gene_compare_enrichment_bysex="mixed";
else if in1 and in3 then   gene_compare_enrichment_bysex="male_both";
else if in2 and in3 then  gene_compare_enrichment_bysex="female_both";
else if in3 then gene_compare_enrichment_bysex="all_both";
else if in1 then gene_compare_enrichment_bysex="male_only";
else  gene_compare_enrichment_bysex="female_only";
drop count percent;
run;

proc freq data=gene_enrich_by_sex;
tables gene_compare_enrichment_bysex;
run;

/*
all_both 287 
female_both 10 
female_only 123 
male_both 29 
male_female_both 1 
male_only 575 
mixed 12 
*/

/*no mutigene
all_both 245 
female_both 10 
female_only 112  
male_both 29 
male_female_both 1  
male_only 476 
mixed 11 
*/




