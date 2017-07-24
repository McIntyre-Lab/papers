ods listing; ods html close;
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
libname event '!MCLAB/event_analysis/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname conesa '!MCLAB/conesa_pacbio/sas_data';
libname refseq '!MCLAB/event_analysis/refseq_fusions/sas_data';

/* Create a list of transcripts from the genes that are expressed */


data exp_genes;
  set event.flag_gene_expressed;
  drop num_fusion_on;
run;

data xs2gene;
  set event.feature2xs2gene_nomulti;
  drop feature_id;
run;

proc sort data=exp_genes;
  by gene_id;
proc sort data=xs2gene nodup;
  by gene_id transcript_id;
run;

data xs2gene_on;
  merge exp_genes (in=in1) xs2gene (in=in2);
  by gene_id;
  if in1 and in2;
  flag_xscript_gene_exp=flag_gene_expressed;
  keep transcript_id flag_xscript_gene_exp flag_gene_only_mult;
run;

data xscript_list;
  set xs2gene;
  keep transcript_id;
run;

proc sort data=xs2gene_on;
  by transcript_id;
proc sort data=xscript_list;
  by transcript_id;
run;

data flag_xscript_w_gene_on;
   merge xscript_list (in=in1) xs2gene_on (in=in2);
   by transcript_id;
   if not in1 then flag_xscript_gene_exp=0;
run;

proc freq data=flag_xscript_w_gene_on;
  tables flag_xscript_gene_exp*flag_gene_only_mult;
run;

/*

flag_xscript_gene_exp
          flag_gene_only_mult

Frequency|
Percent  |
Row Pct  |
Col Pct  |       0|  Total
---------+--------+
       0 |  15892 |  15892
         |  17.31 |  17.31
         | 100.00 |
         |  17.31 |
---------+--------+
       1 |  75926 |  75926
         |  82.69 |  82.69
         | 100.00 |
         |  82.69 |
---------+--------+
Total       91818    91818
           100.00   100.00


75926 of 91818 transcripts where gene is expressed
*/


/* Make permenant */

data event.flag_xscript_w_gene_on;
   set flag_xscript_w_gene_on;
run;

