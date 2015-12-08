/* Importing transcript-annotated junctions into SAS and assemble database */

libname splice '/home/jrbnewman/McLab/junction_annotations/sas_data/';

    data WORK.XSCRIPT_JUNCTIONS    ;
    %let _EFIERR_ = 0; /* set the ERROR detection macro variable */
    infile '/home/jrbnewman/McLab/junction_annotations/generated_files/hg19_xscript_junctions.csv' delimiter = ',' MISSOVER DSD lrecl=32767 firstobs=2 ;
       informat junction_id $50. ;
       informat xscript $50. ;
       informat gene $50. ;
       format junction_id $50. ;
       format xscript $50. ;
       format gene $50. ;
    input
                junction_id $
                xscript $
                gene $
    ;
    if _ERROR_ then call symputx('_EFIERR_',1);  /* set ERROR detection macro variable */
    run;

/* want to cat all xscripts and genes by junction_id so need to separate these datasets */

data junction_gene;
   set xscript_junctions;
   keep junction_id gene;
run;

data junction_xscript;
   set xscript_junctions;
   keep junction_id xscript;
run;


/*sort and remove duplicate observations */

proc sort data=junction_xscript nodup;
   by junction_id xscript;
run; *0 removed;

proc sort data=junction_gene nodup;
   by junction_id gene;
run; *378551 dups removed, 60617 remaining;

/* get counts */

proc freq noprint data=junction_xscript;
   tables junction_id / out=junction_xscript_cnt;
run;

proc sort data=junction_xscript_cnt;
  by descending count;
run;
*max=58 xscripts per junction;


proc freq noprint data=junction_gene;
   tables junction_id / out=junction_gene_cnt;
run;

proc sort data=junction_gene_cnt;
  by descending count;
run;
*max=2 genes per junction; *flag multigene later!;


/* cat xscripts or genes into a single cell per junctions */
/* do for junctions by xscript */

data junction_xscript_cat; 
  array xscripts[58] $ 50;

  retain xscripts1-xscripts58;

  set junction_xscript;
  by junction_id;
  
  if first.junction_id then do;
     call missing(of xscripts1-xscripts58);
     records = 0;
  end;

  records + 1;
  xscripts[records]=xscript;
  if last.junction_id then output;
run;

  *clean up the output file;

data junction_xscript_cat2;
  set junction_xscript_cat;
  length xscripts_cat $ 3000;
  rename records= num_xscripts;
         xscripts_cat= catx("|", OF xscripts1-xscripts58);
  drop xscripts1-xscripts58 xscript ;
  run;

*199181 obs remaining;

/* now do for junctions by gene */

data junction_gene_cat; 
  array genes[2] $ 50;

  retain genes1-genes2;

  set junction_gene;
  by junction_id;
  
  if first.junction_id then do;
     call missing(of genes1-genes2);
     records = 0;
  end;

  records + 1;
  genes[records]=gene;
  if last.junction_id then output;
run;

  *clean up the output file;

data junction_gene_cat2;
  set junction_gene_cat;
  length genes_cat $ 120;
  rename records= num_genes;
         genes_cat= catx("|", OF genes1-genes2);
  drop genes1-genes2 gene ;
  run;
*199181 obs remaining;

/* merge to make database of transcript-annotated junctions */

proc sort data=junction_xscript_cat2;
   by junction_id;
run;

proc sort data=junction_gene_cat2;
   by junction_id;
run;

data junctions_annotated oops_xscript oops_gene;
   merge junction_xscript_cat2 (in=in1) junction_gene_cat2 (in=in2);
   by junction_id;
   if in1 and in2 then output junctions_annotated;
   else if in1 then output oops_xscript; *0 obs, yay;
   else output oops_gene; *0 obs, yay;
run;

/* make data permenant */

data splice.junctions_xscript_annotated_hg19;
   set junctions_annotated;
   if num_genes>1 then flag_multigene=1; *flag for multigene;
   else flag_multigene=0;
run;


/* Clean up */
    proc datasets nolist;
        delete xscript_junctions junction_gene junction_xscript
        junction_xscript_cnt junction_gene_cnt junction_xscript_cat
        junction_xscript_cat2 junction_gene_cat junction_gene_cat2
        junctions_annotated oops_xscript oops_gene
        ;
        run;
        quit;
