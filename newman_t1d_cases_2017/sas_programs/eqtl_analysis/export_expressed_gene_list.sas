Update git
	- check documentation is up to data
	- add an "export expressed immunogenes" step for eQTLs
	- what files do I need?
		- 

/* set libraries */

libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname splicing '/mnt/data/splicing';

/* Get list of genes with exon detection */

data exp_fusions;
   set con.fusions_on_gt_apn0;
   if flag_cd19_on=1 or flag_cd8_on=1 or flag_cd4_on=1;
   keep fusion_id;
run;

data fus2gene;
    set fus.hg19_aceview_fusions_si_info;
    keep fusion_id gene_id;
run;

proc sort data=exp_fusions;
   by fusion_id;
proc sort data=fus2gene nodup;
   by fusion_id gene_id;
run;

data exp_fus2gene;
   merge exp_fusions (in=in1) fus2gene (in=in2);
   by fusion_id;
   if in1 and in2 then output;
run;

data gene_exp_fusion;
   set exp_fus2gene;
   keep gene_id;
run;

proc sort data=gene_exp_fusion nodup;
   by gene_id;
run;


/* Get list of genes with splicing detection */

data gene_exp_splicing;
   set splicing.splicing_results_clean;
   keep gene_id;
run;

proc sort data=gene_exp_splicing nodup;
   by gene_id;
run;

/* Stack lists, drop duplicates */

data genes_exp;
   set gene_exp_splicing gene_exp_fusion;
run;

proc sort data=genes_exp nodup;
    by gene_id;
run;

/* Extract autoimmune genes only */

data immunogenes;
   set con.immunogene_flags;
   if flag_immuno_gene=1;
   keep gene_id;
run;


proc sort data=immunogenes;
   by gene_id;
proc sort data=genes_exp;
   by gene_id;
run;

data immunogenes_w_exp no_flags no_exp;
   merge genes_exp (in=in1) immunogenes (in=in2);
   by gene_id;
   if in1 and in2 then output immunogenes_w_exp;
   else if in1 then output no_flags;
   else output no_exp;
run;

/* Export list of autoimmune genes for pulling SNPs */

proc export data=immunogenes_w_exp outfile='/home/jrbnewman/concannon/eqtl_analysis/immunogenes.txt' dbms=tab replace; putnames=no; run;


