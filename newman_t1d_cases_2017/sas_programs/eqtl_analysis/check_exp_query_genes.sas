/* Check fusion expression of (if no fusion, then NO SPLICING!):
      TYK2, IFIH1, IL2RA, GLIS3, INS-IGF2, BCAR1, GPR183

 Also check SNP selection

*/

libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';
libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

data expressed_fusions;
   set con.fusions_on_gt_apn0;
   if flag_cd19_on=1 or flag_cd4_on=1 or flag_cd8_on=1;
   keep fusion_id;
run;

data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   if gene_id='TYK2'	SNPs, no LD?
   or gene_id='IFIH1'	SNPs, no LD?
   or gene_id='IL2RA'	SNPs, no LD?
   or gene_id='GLIS3'	No data!
   or gene_id='GPR183'	NO SNPS
   or gene_id='INS-IGF2'	No data!
   or gene_id='BCAR1';	No data!
   keep fusion_id gene_id;
run;

proc sort data=expressed_fusions;
   by fusion_id;
proc sort data=fus2gene nodup;
   by fusion_id;
run;

data query_genes_w_exp no_exp;
   merge fus2gene (in=in1) expressed_fusions (in=in2);
   by fusion_id;
   if in1 and in2 then output query_genes_w_exp;
   else if in1 then output no_exp;
run;


data query_genes_w_exp2;
   set query_genes_w_exp;
   keep gene_id;
run;

data query_genes_w_no_exp;
   set no_exp;
   keep gene_id;
run;

proc sort data=query_genes_w_exp2 nodup;
   by gene_id;
proc sort data=query_genes_w_no_exp nodup;
   by gene_id;
run;

data query_genes_exp query_genes_nonexp;
   merge query_genes_w_exp2 (in=in1) query_genes_w_no_exp (in=in2);
   by gene_id;
   if in1 then output query_genes_exp;
   else output query_genes_nonexp;
run;

* all query genes have expression. Do they have SNPs?;

/

