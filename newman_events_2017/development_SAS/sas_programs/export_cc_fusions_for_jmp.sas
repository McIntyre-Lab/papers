libname cc '!PATCON/case_control/sas_data';
libname con '!PATCON/sas_data';
libname cclocal '/mnt/store/cc_sandbox/sas_data';
libname hg19 '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname sem '!MCLAB/grants/ase_grn/INSR-mTOR_T1D/sas_data';

/* Export counts for all immunogene fusions so I can try to cluster samples to ID cell types
For both JMP and Galaxy, I need a wide dataset and a design file
*/

data immunogenes;
  set con.immunobase_gene_flags;
  if flag_immuno_gene=1 then output;
  if gene_id="CD4" then output;
  if gene_id="CD19" then output;
  if gene_id="CD27" then output;
  if gene_id="CD25" then output;
  if gene_id="CD3DandMPZL2" then output;
  if gene_id="CD3EAP" then output;
  if gene_id="CD3E" then output;
  if gene_id="CD3G" then output;
  if gene_id="FOXP3" then output;
  if gene_id="IKZF1" then output;
  if gene_id="IKZF2" then output;
  if gene_id="IKZF3" then output;
  if gene_id="IKZF4" then output;
  if gene_id="IKZF5" then output;
  keep gene_id;
run;

data fus2gene;
  set hg19.hg19_aceview_fusions_si_info;
  where flag_multigene=0;
  keep fusion_id gene_id ;
run;

proc sort data=immunogenes nodup;
  by gene_id;
proc sort data=fus2gene nodup;
  by gene_id fusion_id;
run;

data immuno_fus;
  merge immunogenes (in=in1) fus2gene (in=in2);
  by gene_id;
  if in1 and in2;
run;

data counts;
  set cc.cc_counts_w_zeros;
  name2=name+0;
  drop name;
  rename name2=name;
run;

data immuno_fus2;
  set immuno_fus;
  keep fusion_id;
run;

proc sort data=immuno_fus2 nodup;
   by fusion_id;
proc sort data=counts;
   by fusion_id;
run;

data fus_w_counts;
  merge immuno_fus2 (in=in1) fus_w_counts (in=in2);
  by fusion_id;
  if in1 and in2;
run;




