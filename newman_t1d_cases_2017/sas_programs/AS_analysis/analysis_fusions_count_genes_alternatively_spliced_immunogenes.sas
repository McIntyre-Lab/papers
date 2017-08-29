
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


data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=ai nodup;
  by gene_id;
proc sort data=t1d nodup;
  by gene_id;
proc sort data=fusion_spec_w_gene;
  by gene_id;
run;

data fusion_spec_w_gene_immuno;
   merge fusion_spec_w_gene (in=in1) ai (in=in2) t1d (in=in3);;
   by gene_id;
   if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
   if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in1 then output;
run;


/* Count DD exons in genes expressed in both cell types! */


%macro genecounts(geneflag);

%macro count_dd(cell1,cell2);

data &cell1.&cell2._fus;
   set fusion_spec_w_gene_immuno;
   where &geneflag.=1;
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
run; 

%mend;

%genecounts(flag_immuno_gene);
%genecounts(flag_immunobase_diabetes_gene);
/* RESULTS:

AUTOIMMUNE GENES

CD4 v CD8:
      flag_gene_                             Cumulative    Cumulative
       cd4cd8_as    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        1306       78.67          1306        78.67
               1         354       21.33          1660       100.00


CD4 v CD19:
      flag_gene_                             Cumulative    Cumulative
      cd4cd19_as    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0        1105       67.50          1105        67.50
               1         532       32.50          1637       100.00


CD8 v CD19:


       flag_gene_                             Cumulative    Cumulative
       cd8cd19_as    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0        1123       68.39          1123        68.39
                1         519       31.61          1642       100.00


617 genes total AS


T1D GENES

CD4 v CD8:
       flag_gene_                             Cumulative    Cumulative
        cd4cd8_as    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0         330       82.71           330        82.71
                1          69       17.29           399       100.00

                             The SAS System              10:49 Friday, May 5, 20

CD4 v CD19:
    flag_gene_                             Cumulative    Cumulative
    cd4cd19_as    Frequency     Percent     Frequency      Percent
-------------------------------------------------------------------
             0         276       70.59           276        70.59
             1         115       29.41           391       100.00


CD8 v CD19:

  flag_gene_                             Cumulative    Cumulative
  cd8cd19_as    Frequency     Percent     Frequency      Percent
-----------------------------------------------------------------
           0         283       72.01           283        72.01
           1         110       27.99           393       100.00

133 genes total AS

*/





