/* Create gene-level datasets for MMC (average across on fusions)
   I am also going to iterate through various options in my MMC gene list */

ods listing; ods html close;
libname cc '!PATCON/case_control/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

/* options for genelist:
ds487		487 genes that are highly diff spliced between cases and controls, CD4+/CD25- and CD4+/CD25+
ds487de		487 genes plus genes that are DE for exons
ds487ds		487 genes plus genes that are differentially spliced
ds487deds	487 genes plus genes that are differentially spliced and exons DE
deds		Genes that are differentially spliced and exons DE
de		Genes that are DE
ds		Genes that are differentially spliced

For each, I will make a control-only, case-only and together dataset.
  For DE lists, DE in CD4+/CD25+, CD4+/CD25- and DE in both
*/


data gene_list;
  set cc.gene_list_for_mmc;
run;

/* 487 genes + DE genes + DS genes */

    /* DE/DD exons AND DE/DD events */



*CD425L DE genes;

data gene_list2;
   set gene_list;
   if flag_ds_425v426=1 or (flag_gene_fus_de_cd425l=1 and flag_gene_as_de_cd425l=1);
   keep gene_id;
run;

%makemmc(487_ds_genes_cd425l_de_genes_fus_and_as,all);
%makemmc(487_ds_genes_cd425l_de_genes_fus_and_as,controls);
%makemmc(487_ds_genes_cd425l_de_genes_fus_and_as,cases);

*CD425L DD genes;

data gene_list2;
   set gene_list;
   if flag_ds_425v426=1 or (flag_gene_fus_dd_cd425l=1 and flag_gene_as_dd_cd425l=1);
   keep gene_id;
run;

%makemmc(487_ds_genes_cd425l_dd_genes_fus_and_as,all);
%makemmc(487_ds_genes_cd425l_dd_genes_fus_and_as,controls);
%makemmc(487_ds_genes_cd425l_dd_genes_fus_and_as,cases);

*CD425H DE genes;

data gene_list2;
   set gene_list;
   if flag_ds_425v426=1 or (flag_gene_fus_de_cd425h=1 and flag_gene_as_de_cd425h=1);
   keep gene_id;
run;

%makemmc(487_ds_genes_cd425h_de_genes_fus_and_as,all);
%makemmc(487_ds_genes_cd425h_de_genes_fus_and_as,controls);
%makemmc(487_ds_genes_cd425h_de_genes_fus_and_as,cases);

*CD425H DD genes;

data gene_list2;
   set gene_list;
   if flag_ds_425v426=1 or (flag_gene_fus_dd_cd425h=1  and flag_gene_as_dd_cd425h=1);
   keep gene_id;
run;

%makemmc(487_ds_genes_cd425h_dd_genes_fus_and_as,all);
%makemmc(487_ds_genes_cd425h_dd_genes_fus_and_as,controls);
%makemmc(487_ds_genes_cd425h_dd_genes_fus_and_as,cases);

*CD425L+CD425H DE genes;

data gene_list2;
   set gene_list;
   if flag_ds_425v426=1 or (flag_gene_fus_de_cd425l=1 and flag_gene_fus_de_cd425h=1
    and flag_gene_as_de_cd425l=1 and flag_gene_as_de_cd425h=1);
   keep gene_id;
run;

%makemmc(487_ds_genes_cd425l_cd425h_de_genes_fus_and_as,all);
%makemmc(487_ds_genes_cd425l_cd425h_de_genes_fus_and_as,controls);
%makemmc(487_ds_genes_cd425l_cd425h_de_genes_fus_and_as,cases);


*CD425L+CD425H DD genes;

data gene_list2;
   set gene_list;
   if flag_ds_425v426=1 or (flag_gene_fus_dd_cd425l=1 and flag_gene_fus_dd_cd425h=1
    and flag_gene_as_dd_cd425l=1 and flag_gene_as_dd_cd425h=1) ;
   keep gene_id;
run;

%makemmc(487_ds_genes_cd425l_cd425h_dd_genes_fus_and_as,all);
%makemmc(487_ds_genes_cd425l_cd425h_dd_genes_fus_and_as,controls);
%makemmc(487_ds_genes_cd425l_cd425h_dd_genes_fus_and_as,cases);

