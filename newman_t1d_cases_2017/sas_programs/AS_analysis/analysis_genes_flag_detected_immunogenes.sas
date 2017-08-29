ods listing; ods html close;
libname con '!PATCON/sas_data';
libname fus '!PATCON/useful_human_data/aceview_hg19/fusions/sas_data';
libname splicing '/mnt/data/splicing';

/* Count T1D and AI genes expressed by cell type */



data flag_gene_on;
   set con.flag_gene_detection_by_cell;
run;

data ai t1d;
   set con.immunobase_gene_flags;
   if flag_immuno_gene=1 then output ai;
   if flag_immunobase_diabetes_gene=1 then output t1d;
   keep gene_id;
run;

proc sort data=flag_gene_on nodup;
   by gene_id;
proc sort data=ai nodup;
   by gene_id;
proc sort data=t1d nodup;
   by gene_id;
run;

data genes_on_immuno;
   merge flag_gene_on (in=in1) ai (in=in2) t1d (in=in3);
   by gene_id;
   if in2 then flag_immuno_gene=1; else flag_immuno_gene=0;
   if in3 then flag_immunobase_diabetes_gene=1; else flag_immunobase_diabetes_gene=0;
   if in1 then output;
run;

/* Count autoimmune genes on */

proc freq data=genes_on_immuno noprint;
   where flag_immuno_gene=1;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on / out=ai_gene_count;
proc print data=ai_gene_count;
run;


/*
                            flag_
 flag_cd4_    flag_cd8_     cd19_
  gene_on      gene_on     gene_on    COUNT

     0            0           0         162
     0            0           1          16
     0            1           0           6
     0            1           1           6
     1            0           0           7
     1            0           1           2
     1            1           0          25
     1            1           1        1628

*/

/* Count T1D genes on */


proc freq data=genes_on_immuno noprint;
   where flag_immunobase_diabetes_gene=1;
   tables flag_cd4_gene_on*flag_cd8_gene_on*flag_cd19_gene_on / out=t1d_gene_count;
proc print data=t1d_gene_count;
run;

/*
                            flag_
 flag_cd4_    flag_cd8_     cd19_
  gene_on      gene_on     gene_on    COUNT

     0            0           0         27
     0            0           1          5
     0            1           1          1
     1            0           0          2
     1            1           0          8
     1            1           1        389

*/

proc export data=ai_gene_count
     outfile="!PATCON/pipeline_output/differential_splicing_counts/autoimmune_gene_counts.csv"
     dbms=csv replace;
run;

proc export data=t1d_gene_count
     outfile="!PATCON/pipeline_output/differential_splicing_counts/diabetes_gene_counts.csv"
     dbms=csv replace;
run;


