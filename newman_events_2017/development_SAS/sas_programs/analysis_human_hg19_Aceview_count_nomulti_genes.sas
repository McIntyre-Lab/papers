ods listing; ods html close;
libname con '!PATCON/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';

* genes on in all cell types;


data genes_on;
  set con.flag_gene_detection_by_cell;
  where flag_cd4_gene_on=1 and flag_cd8_gene_on=1 and flag_cd19_gene_on=1;
  keep gene_id;
run;

data genes_multi;
   set hg19.hg19_aceview_fusions_si_info;
   where flag_multigene=1;
   keep gene_id;
run;

proc sort data=genes_on;
   by gene_id;
proc sort data=genes_multi nodup; 
  by gene_id;
run;

22670

data multi nomulti;
  merge genes_on (in=in1) genes_multi (in=in2);
  by gene_id;
  if in1 and in2 then output multi;
  else if in1 then output nomulti;
run;

/*
17791 with multi
23236 without multi
*/
