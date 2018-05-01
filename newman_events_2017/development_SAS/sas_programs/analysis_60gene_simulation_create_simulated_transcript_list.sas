/* Libraries */
libname event  '!MCLAB/event_analysis/sas_data';
libname eventloc '/mnt/store/event_sandbox/sas_data';
libname evspl '!MCLAB/conesa_pacbio/sas_data/splicing';
libname mm10 '!MCLAB/useful_mouse_data/mm10/sas_data';
ods listing; ods html close;

/* 60 gene simulation, create a gene and transcript list that I can refer back to */

data xs2gene;
  set event.feature2xs2gene;
  keep transcript_id gene_id;
run;

data genelist_1;
  format gene_id $11.;
  input gene_id $;
  datalines;
  101214
  106581
  108011
  114565
  114863
  11864
  12419
  12830
  142980
  16150
  16468
  170787
  17350
  20509
  217378
  22029
  22276
  225608
  22654
  228960
  231386
  232337
  234366
  234725
  269120
  27416
  28240
  56812
  58175
  66433
  66855
  66940
  67669
  69930
  70380
  70382
  70796
  71893
  72205
  72344
  72999
  73750
  74080
  74140
  74190
  76577
  77613
  77733
  94281
  ;
run;

data genelist_2;
  set event.genes_w_nic_junction_10genes;
run;

proc sort data=genelist_1;
  by gene_id;
proc sort data=genelist_2;
  by gene_id;
run;

data genelist_for_sim;
  merge genelist_1 (in=in1) genelist_2 (in=in2);
  by gene_id;
run;

proc sort data=genelist_for_sim;
  by gene_id;
proc sort data=xs2gene nodup;
  by gene_id;
run;

data xslist_sim;
  merge genelist_for_sim (in=in1) xs2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

/* MAke permenant */

data event.polyester_xs_list_60genes;
  set xslist_sim;
run;

