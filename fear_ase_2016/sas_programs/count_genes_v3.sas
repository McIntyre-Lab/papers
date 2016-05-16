libname dmel "McLab/useful_dmel_data/flybase551/sasdata";
libname cegs "McLab/cegs_ase_paper/sas_data";

libname dmel "/home/mcintyre/S/SHARE/McIntyre_Lab/useful_dmel_data/flybase551/sasdata";
libname cegs "/home/mcintyre/S/SHARE/McIntyre_Lab/cegs_ase_paper/sas_data";

proc contents data=cegs.Cis_est_v13;
run;

proc sort data=cegs.Cis_est_v13;
by fusion_id;
run;

data sig_est;
set cegs.Cis_est_v13;
if flag_all_ai=1;
run;

proc sort data=dmel.Fb551_si_fusions_unique_flagged;
by fusion_id;

*count direction of AI within a line mating status;
data sig_anno;
merge sig_est (in=in1) dmel.Fb551_si_fusions_unique_flagged(in=in2);
by fusion_id;
if in1;
if q5_mean_theta>0.5 then direction=1;
else direction=-1;
run;


proc sort data= sig_anno;
by line mating_status;

proc freq data=sig_anno noprint;
by line mating_status;
tables symbol_cat*direction/out=count_direction;
run;


proc freq data=count_direction noprint;
by line mating_status;
tables symbol_cat/out=count_switch;
run;



proc freq data=count_switch;
tables count;
run;

*perhaps we have ~3% ase specific splicing ;
*need to make wiggles;

data check_switch;
set count_switch;
if count =2;
keep symbol_cat line mating_status;
run;

data cegs.ai_switch;
set check_switch;
run;

*how many unqiue genes are potential splicers;

proc freq data=check_switch noprint;
tables symbol_cat/out=genes_with_splice_ai;
run;

data cegs.genes_qith_spice_ai;
set genes_with_splice_ai;
rename count=num_lines_env;
run;
proc sort data=check_switch;
by symbol_cat ;
*line mating_status;
run;

proc sort data=sig_anno;
by symbol_cat line mating_status;

data find_switch;
merge sig_anno (in=in1) check_switch(in=in2);
by symbol_cat line mating_status;
if in2 then flag_splice_Ai=1;
else flag_splice_ai=0;
run;

*look at ald!;

*should set to non_sig and look at lists;

*how many genes;
proc freq data=find_switch noprint;
tables symbol_cat/out=count_genes_with_ai;
run;

*2252 genes with AI;
proc freq data=count_genes_with_ai;
tables count;
run;
*12.57 % in only one line_mating_status_exon;
*88% of genes have at least 2 ;
*76% have at least 3;

data lots_of_ai;
set count_genes_with_ai;
if count ge 50;
keep symbol_cat;


proc sort data=lots_of_ai;
by symbol_cat ;
run;

proc sort data=find_switch;
by symbol_cat ;

data genes_lots_ai;
merge lots_of_ai (in=in1) find_switch(in=in2);
by symbol_cat ;
if in1 then flag_lots_Ai=1;
else flag_lots_ai=0;
run;

proc freq data=genes_lots_ai;
tables flag_lots_ai*flag_splice_ai;
run;

data peek;
set genes_lots_ai ;
where  flag_splice_ai=1;
run;

proc sort data=genes_lots_ai;
by symbol_cat mating_status;

proc freq data=genes_lots_ai noprint;
by symbol_cat mating_status;
tables line*direction /out=count_line_dir;
run;

proc freq data=count_line_dir noprint;
by symbol_cat mating_status;
tables direction/out=count_dir;
run;

proc sort data=count_dir;
by mating_status;
run;

proc freq data=count_dir noprint;
by  mating_status;
tables symbol_cat/out=count_gene_ai_direction;
run;

data flag_gene_ai_direction;
set count_gene_ai_direction;
if count=1 then gene_direction="single";
else if count=2 then gene_direction="mixed";
keep symbol_cat mating_status gene_direction;
run;

proc freq data=flag_gene_ai_direction;
tables gene_direction;
run;

proc sort data=genes_lots_ai;
by symbol_cat mating_status;

proc sort data=flag_gene_ai_direction;
by symbol_cat mating_status;

data genes_with_ai;
merge genes_lots_ai (in=in1) flag_gene_ai_direction (in=in2);
by symbol_cat mating_status;
run;
