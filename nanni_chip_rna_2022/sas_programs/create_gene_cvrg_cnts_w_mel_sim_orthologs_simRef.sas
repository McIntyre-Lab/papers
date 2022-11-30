libname ortho "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/sas_data/RNAseq_ortho";
libname df  "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/svn_lmm_dros_head_data/ethanol_srna/sas_data";

libname drosRNA "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/sas_data/RNAseq";

/*
combine mel and sim gene level coverage counts using ortholog file


Ortholog file (mel_geneID, sim_geneID): /blue/mcintyre/share/etoh_srna/ortholog_files/dmel_orthologs_dsim_fb_2017_04.csv
    scp above file to /nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq


*/



/* create species stacks */
data cc_mel_stack ;
set ortho.cvrg_gene_wt_log_uq_apn_simRef ;
if find(sampleID, "mel") ge 1 ;
rename geneID = sim_geneID ; /* mel aligned to sim ref */
run ;

data cc_sim_stack ;
set ortho.cvrg_gene_wt_log_uq_apn_simRef ;
if find(sampleID, "sim") ge 1 ;
rename geneID = sim_geneID ; /* sim aligned to sim ref */
run ;

proc freq data = cc_mel_stack ;
tables sim_geneID / out = cnts_mel ;
run; /* 11261 genes, 24 obs per gene  */
proc freq data = cc_sim_stack ;
tables  sim_geneID / out = cnts_sim ;
run; /* 11261 genes, 23obs per gene (dropped 1 sample */


/* import ortho table mel-2-sim */
proc import datafile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq/dmel_orthologs_dsim_fb_2017_04.csv"
out = ortho2
dbms = csv replace ;
guessingrows = MAX ;
run;

proc contents data = ortho2 ; run;

proc freq data = ortho2 ;
tables mel_geneID / out = ortho_mel_cnts;
run;
data ck_mel ;
set ortho_mel_cnts ;
where count ne 1 ;
run ;
    /* 14,601 obs, 13,373 mel genes, 779 where count > 1 */

proc freq data = ortho2 ;
tables sim_geneID / out = ortho_sim_cnts ;
run;

data ck_sim ;
set ortho_sim_cnts ;
where count ne 1 ;
run ;
    /* 14,601 obs, 13,306 sim genes, 724 where count > 1 */

data ortho ;
set ortho2 ;
where flag_one2one_ortholog = 1 ;
melID_simID = compress(mel_geneID||'_'||sim_geneID) ;
keep melID_simID mel_geneID sim_geneID;
run ; 

proc sort data = ortho nodups ;
by _all_ ;
run; /* 12386 obs */

/* sort and merge ortho to cc mel stack */
proc sort data = ortho ;
by sim_geneID ;
proc sort data = cc_mel_stack ;
by sim_geneID;
run;

data cc_mel_ortho ;
merge cc_mel_stack (in=in1) ortho (in=in2) ;
by sim_geneID ;
if in1 ;
run;

data cc2_mel_ortho ;
retain melID_simID mel_geneID sim_geneID  ;
set cc_mel_ortho ;
drop mel_geneID sim_geneID ;
where melID_simID ne "" ;
run;

proc freq data = cc2_mel_ortho ;
table melID_simID / out = cnt_mel_mel_sim;
run ;
    /* 9806 melID_simID */

/* sort and merge ortho to cc sim stack */
proc sort data = ortho ;
by sim_geneID ;
proc sort data = cc_sim_stack ;
by sim_geneID ;
run ;

data cc_sim_ortho ;
merge cc_sim_stack (in=in1) ortho (in=in2) ;
if in1 ;
by sim_geneID ;
run;

data cc2_sim_ortho ;
retain melID_simID mel_geneID sim_geneID  ;
set cc_sim_ortho ;
drop mel_geneID sim_geneID ;
where melID_simID ne "" ;
run;

proc freq data = cc2_sim_ortho ;
table melID_simID / out = cnt_sim_mel_sim;
run ;
    /* 9806 melID_simID */

data cvrg_gene_mel_sim_orthologs ;
set cc2_mel_ortho cc2_sim_ortho ;
run ;

proc freq data = cvrg_gene_mel_sim_orthologs ;
tables melID_simID / out = cnt_ID  ;
run;  /* 47 samples per ID */


data ortho.cvrg_norm_gene_ortho_2_sim ;  /* stacked dataset */
set cvrg_gene_mel_sim_orthologs ;
run;


/* make wide for log_ug_apn and log_ug_apn_noZero */
proc contents data = ortho.cvrg_norm_gene_ortho_2_sim ; run;

proc sort data = ortho.cvrg_norm_gene_ortho_2_sim ; 
by melID_simID ;
run;

%macro flip (var, shrt) ;

proc transpose data =ortho.cvrg_norm_gene_ortho_2_sim out = wide_&var. ; 
by melID_simID ;
id sampleID ;
var &var. ;
run;

data ortho.cvrGene_ortho2sim_&shrt. ;
set wide_&var. ;
drop _name_ ;
run ;

proc export data = ortho.cvrGene_ortho2sim_&shrt.
outfile = "/nfshome/ammorse/mclab/SHARE/McIntyre_Lab/Dros_CHIP_RNA_ms/RNAseq/Coverage_counts/gene_cvrg_cnts_ortho/cvrGene_ortho2sim_&shrt..csv"
dbms = csv replace;
run ;

%mend ;

%flip (log_uq_apn, log_uq_apn) ;
%flip (log_uq_apn_noZero, log_uq_apn_no0) ;
 

