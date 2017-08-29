
/* Libraries */

libname eqtl '/home/jrbnewman/concannon/eqtl_analysis/sas_data';
libname splicing '/mnt/data/splicing';
libname splice '/mnt/data/splice';
libname con '/home/jrbnewman/concannon/sas_data';
libname fus '/home/jrbnewman/McLab/useful_human_data/aceview_hg19/fusions/sas_data';

/* I want to now count the number of genes detected in each cell type pair that have at least 1 fusion
   differentially detected. This will be evidence of alternative splice (gene must be detected in both cell types!)
*/

data flag_splicing_spec ;
   set splicing.flag_splicing_by_gene_dtct;
run;

data event2gene;
  set splice.splicing_events_annotations;
  keep event_id gene_id;
run;


proc sort data=flag_splicing_spec;
   by event_id;
proc sort data=event2gene nodup;
   by event_id gene_id;
run;

data splicing_spec_w_gene;
  merge flag_splicing_spec (in=in1) event2gene (in=in2);
  by event_id;
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
proc sort data=splicing_spec_w_gene;
  by gene_id;
run;

data splicing_spec_w_gene_immuno;
  merge splicing_spec_w_gene (in=in1) ai (in=in2) t1d (in=in3);
  by gene_id;
  if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
  if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
  if in1 then output;
run;


/* Count DD exons in genes expressed in both cell types! */


%macro genecounts(geneflag);

%macro count_dd(cell1,cell2);

data &cell1.&cell2._fus;
   set splicing_spec_w_gene_immuno;
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
   output out=&cell1.&cell2._num_dd sum=num_events_&cell1.&cell2._dd;
run;

data flag_gene_dd_&cell1.&cell2.;
   set &cell1.&cell2._num_dd;
   if num_events_&cell1.&cell2._dd > 0 then flag_gene_&cell1.&cell2._as=1;
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

data cd4cd8_as;
   set flag_gene_dd_cd4cd8;
   where flag_gene_cd4cd8_as=1;
   keep gene_id;
run;

data cd4cd19_as;
   set flag_gene_dd_cd4cd19;
   where flag_gene_cd4cd19_as=1;
   keep gene_id;
run;

data cd8cd19_as;
   set flag_gene_dd_cd8cd19;
   where flag_gene_cd8cd19_as=1;
   keep gene_id;
run;

data genes_as;
  set cd4cd8_as cd4cd19_as cd8cd19_as;
run;

proc sort data=genes_as nodup;
  by gene_id;
run; *11627 genes;


%mend;

%genecounts(flag_immuno_gene);
%genecounts(flag_immunobase_diabetes_gene);

/* RESULTS:

AUTOIMMUNE:
CD4 v CD8:
        flag_gene_                             Cumulative    Cumulative
         cd4cd8_as    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0         678       46.06           678        46.06
                 1         794       53.94          1472       100.00

CD4 v CD19:
     flag_gene_                             Cumulative    Cumulative
     cd4cd19_as    Frequency     Percent     Frequency      Percent
--------------------------------------------------------------------
              0         457       31.11           457        31.11
              1        1012       68.89          1469       100.00

CD8 v CD19:

          flag_gene_                             Cumulative    Cumulative
          cd8cd19_as    Frequency     Percent     Frequency      Percent
    ---------------------------------------------------------------------
                   0         452       30.73           452        30.73
                   1        1019       69.27          1471       100.00

1111 genes DS

1288 genes in total

T1D:
CD4 v CD8:
      flag_gene_                             Cumulative    Cumulative
       cd4cd8_as    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0         165       47.41           165        47.41
               1         183       52.59           348       100.00

CD4 v CD19:
         flag_gene_                             Cumulative    Cumulative
         cd4cd19_as    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0         108       31.21           108        31.21
                  1         238       68.79           346       100.00

CD8 v CD19:
      flag_gene_                             Cumulative    Cumulative
      cd8cd19_as    Frequency     Percent     Frequency      Percent
---------------------------------------------------------------------
               0         106       30.64           106        30.64
               1         240       69.36           346       100.00

262 genes DS
296 genes total
*/



data genes_as_all;
   set genes_as genes_as_fus;
run;

proc sort data=genes_as_all nodup;
   by gene_id;
run; *15428 genes;

