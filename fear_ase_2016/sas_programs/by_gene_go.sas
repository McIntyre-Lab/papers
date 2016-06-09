

data all;
set  cegs.Clean_ase_stack;
run;

proc sort data=all;
by fusion_id;

proc sort data=fusion2go;
by fusion_id;

data all_anno;
merge all (in=in1) fusions2go(in=in2);
by fusion_id;
if in1;
run;

proc freq data=all_anno noprint;
tables symbol_cat*mating_status/out=genes_tested;
run;

data genes_tested2;
set genes_tested;
rename count=num_lines;
run;

proc sort data=genes_tested2;
by symbol_cat mating_status;
proc sort data=flag_gene_ai_direction ;
by symbol_cat mating_status;

data genes_tested3;
merge genes_tested2 flag_gene_ai_direction (in=in2);
by symbol_cat mating_status;
drop count;
if in2 then flag_sig=1;
else flag_sig=0;
run;

proc sort data=genes_tested3;
by symbol_cat;

proc sort data=genes2go;
by symbol;
*problem;
*ignoring for now;

data check;
set genes_tested3;
where symbol_cat ? "|";
run;

data genes_tested4;
set genes_tested3;
rename symbol_cat=symbol;
run;
data genes_anno;
merge genes_tested4 (in=in1) genes2go;
by symbol;
if in1;
run;

proc sort data=defend2;
by symbol;

data all_anno_defend;
merge genes_anno(in=in1) defend2 ;
by symbol;
if in1;
run;

proc freq data=all_anno_defend;
tables flag_sig*go_defense_551;
run;

