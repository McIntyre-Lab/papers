*libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';

libname eqtl '/mnt/data/eqtls/sas_data';
libname con '/home/jrbnewman/concannon/sas_data';
libname splicing '/mnt/data/splicing/';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/**** Prepare fusion info *****/


/* First get a list of "on" fusions - fusion needs to be on in only one cell type */
data fusions_on;
   set con.fusions_on_gt_apn0;
   if flag_CD19_on=1 or flag_CD4_on=1 or flag_CD8_on=1;
   keep fusion_id flag_CD19_on flag_CD4_on flag_CD8_on;
run;

/* Merge in gene info for fusions */

data fusion2gene;
   set fus.unique_info_fusions_si;
   keep fusion_id gene_id;
run;

proc sort data=fusion2gene;
   by fusion_id;
proc sort data=fusions_on;
   by fusion_id;
run;

data fusion2gene_on;
    merge fusion2gene (in=in1) fusions_on (in=in2);
    by fusion_id;
    if in1 and in2 then output;
run;

/* Get fusion counts */

proc sort data=fusion2gene_on;
   by fusion_id;
run;

proc sort data=con.fusion_q3_norm_data_all;
   by fusion_id;
run;

data fusion_counts;
   merge fusion2gene_on (in=in1) con.fusion_q3_norm_data_all (in=in2);
   by fusion_id;
   if in1 and in2;
run;


/* Merge in design file */

data design_file;
   set con.design_file;
   keep subject_id Name cell_type;
run;

proc sort data=design_file nodup;
   by Name;
run;

proc sort data=fusion_counts;
   by Name;
run;

data fusion_counts_w_info;
    merge design_file (in=in1) fusion_counts (in=in2);
    by Name;
    if in1 and in2 then output;
run;

/* Drop "off" fusions per cell type */

data fusion_counts_w_info2;
    set fusion_counts_w_info;
    if cell_type='CD4' and flag_cd4_on=0 then delete;
    if cell_type='CD8' and flag_cd8_on=0 then delete;
    if cell_type='CD19' and flag_cd19_on=0 then delete;
run;


/* Un-cat gene ids since we will be merging SNPs and fusions on gene later
   For  duplicated entries (since there will be some) we will sort and remove dups later */

data fusion_counts_w_info3(rename=new=gene);
length new $35.; 
set fusion_counts_w_info2; 
do i=1 by 1 while(scan(gene_id,i,'|') ^=' '); 
new=scan(gene_id,i,'|'); 
output; 
end; 
run;


data fusion_counts_w_info4;
    set fusion_counts_w_info3;
    drop gene_id i;
    rename gene=gene_id;
run;



/* Make permenant */
data eqtl.fusion_counts_for_eqtls;
    set fusion_counts_w_info4;
run;

