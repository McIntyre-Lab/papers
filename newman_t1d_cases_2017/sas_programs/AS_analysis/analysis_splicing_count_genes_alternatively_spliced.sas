
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
   set splicing.flag_splicing_by_gene_dtct_v2;
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
  if flag_cd4_gene_on=0 then flag_cd4_on=0;
  if flag_cd8_gene_on=0 then flag_cd8_on=0;
  if flag_cd19_gene_on=0 then flag_cd19_on=0;
run;

/* Count DD exons in genes expressed in both cell types! */

%macro count_dd(cell1,cell2);

data &cell1.&cell2._fus;
   set splicing_spec_w_gene;
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
 

/* RESULTS:

CD4 v CD8:

       flag_gene_                             Cumulative    Cumulative
        cd4cd8_as    Frequency     Percent     Frequency      Percent
 ---------------------------------------------------------------------
                0       10746       59.11         10746        59.11
                1        7434       40.89         18180       100.00

CD4 v CD19:


        flag_gene_                             Cumulative    Cumulative
        cd4cd19_as    Frequency     Percent     Frequency      Percent
  ---------------------------------------------------------------------
                 0        7779       43.27          7779        43.27
                 1       10197       56.73         17976       100.00



CD8 v CD19:

         flag_gene_            18,00                 Cumulative    Cumulative
         cd8cd19_as    Frequency     Percent     Frequency      Percent
   ---------------------------------------------------------------------
                  0        7700       42.78          7700        42.78
                  1       10301       57.22         18001       100.00

*/

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




data genes_as_all;
   set genes_as genes_as_fus;
run;

proc sort data=genes_as_all nodup;
   by gene_id;
run; *15428 genes;


proc freq data=con.gene_diff_splicing_summary;
   where flag_cd4_gene_on=1 and flag_cd19_gene_on=1;
   tables flag_cd4cd19_event_dd;
run;

data check;
  set con.gene_diff_splicing_summary;
  where flag_cd4_gene_on=1 and flag_cd8_gene_on=1;
  if flag_cd4cd8_exon_dd=1 or flag_cd4cd8_event_dd=1;
run;



                  flag_cd4cd19_                             Cumulative    Cumulative
                       event_dd    Frequency     Percent     Frequency      Percent
               ---------------------------------------------------------------------
                              0       31327       75.44         31327        75.44
                              1       10197       24.56         41524       100.00



