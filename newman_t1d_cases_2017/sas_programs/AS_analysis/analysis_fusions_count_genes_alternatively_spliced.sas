
/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* I want to now count the number of genes detected in each cell type pair that have at least 1 fusion
   differentially detected. This will be evidence of alternative splice (gene must be detected in both cell types!)
*/

data flag_fusion_spec ;
   set con.flag_fusions_by_gene_detection2;
run;

data fus2gene;
   set fus.hg19_aceview_fusions_si_info;
   keep fusion_id gene_id;
run;

proc sort data=flag_fusion_spec;
   by fusion_id;
proc sort data=fus2gene nodup;
   by fusion_id gene_id;
run;

data fusion_spec_w_gene;
  merge flag_fusion_spec (in=in1) fus2gene (in=in2);
  by fusion_id;
  if in1 and in2;
run;

/* Count DD exons in genes expressed in both cell types! */

%macro count_dd(cell1,cell2);

data &cell1.&cell2._fus;
   set fusion_spec_w_gene;
   if flag_&cell1._gene_on=1 and flag_&cell2._gene_on=1 then output &cell1.&cell2._fus;
run;


data &cell1.&cell2._flag_dd;
   set &cell1.&cell2._fus;
   if flag_&cell1._on=1 and flag_&cell2._on=0 then flag_&cell1.&cell2._dd=1;
   else if flag_&cell1._on=0 and flag_&cell2._on=1 then flag_&cell1.&cell2._dd=1;
   else if flag_&cell1._on=1 and flag_&cell2._on=1 then flag_&cell1.&cell2._dd=0;
   else if flag_&cell1._on=0 and flag_&cell2._on=0 then flag_&cell1.&cell2._dd=.;
run;

proc sort data=&cell1.&cell2._flag_dd;
   by gene_id;
proc means data=&cell1.&cell2._flag_dd noprint;
   by gene_id;
   var flag_&cell1.&cell2._dd;
   output out=&cell1.&cell2._num_dd sum=num_fusions_&cell1.&cell2._dd;
run;

data flag_gene_dd_&cell1.&cell2.;
   set &cell1.&cell2._num_dd;
   if num_fusions_&cell1.&cell2._dd > 0 then flag_gene_&cell1.&cell2._as=1;
   else flag_gene_&cell1.&cell2._as=0;
run;

proc freq data=flag_gene_dd_&cell1.&cell2.;
   tables flag_gene_&cell1.&cell2._as;
run;
%mend;

%count_dd(cd4,cd8);
%count_dd(cd4,cd19);
%count_dd(cd8,cd19);
 

/* RESULTS:

CD4 v CD8:


 flag_gene_                             Cumulative    Cumulative
  cd4cd8_as    Frequency     Percent     Frequency      Percent
----------------------------------------------------------------
          0       38540       88.98         38540        88.98
          1        4774       11.02         43314       100.00

CD4 v CD19:


     flag_gene_                             Cumulative    Cumulative
     cd4cd19_as    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0       34771       83.00         34771        83.00
              1        7123       17.00         41894       100.00

CD8 v CD19:

      flag_gene_                             Cumulative    Cumulative
      cd8cd19_as    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0       34643       82.73         34643        82.73
               1        7231       17.27         41874       100.00


*/

/* Number of genes alternatively spliced */

data cd4cd8_as_fus;
   set flag_gene_dd_cd4cd8;
   where flag_gene_cd4cd8_as=1;
   keep gene_id;
run;

data cd4cd19_as_fus;
   set flag_gene_dd_cd4cd19;
   where flag_gene_cd4cd19_as=1;
   keep gene_id;
run;

data cd8cd19_as_fus;
   set flag_gene_dd_cd8cd19;
   where flag_gene_cd8cd19_as=1;
   keep gene_id;
run;

data genes_as_fus;
  set cd4cd8_as_fus cd4cd19_as_fus cd8cd19_as_fus;
run;

proc sort data=genes_as_fus nodup;
  by gene_id;
run; *8700 genes;




