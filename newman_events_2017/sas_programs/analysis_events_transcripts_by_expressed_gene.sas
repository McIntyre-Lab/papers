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
  set event.feature2xs2gene;
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
  Col Pct  |       0|       1|  Total
  ---------+--------+--------+
         0 |  22591 |   1157 |  23748
           |  17.56 |   0.90 |  18.46
           |  95.13 |   4.87 |
           |  17.72 | 100.00 |
  ---------+--------+--------+
         1 | 104883 |      0 | 104883
           |  81.54 |   0.00 |  81.54
           | 100.00 |   0.00 |
           |  82.28 |   0.00 |
  ---------+--------+--------+
  Total      127474     1157   128631
              99.10     0.90   100.00



104883 of 128631 transcripts where gene is expressed
*/


/* Make permenant */

data event.flag_xscript_w_gene_on;
   set flag_xscript_w_gene_on;
run;

